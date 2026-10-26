
import streamlit as st
import pandas as pd
import requests
import plotly.express as px
import plotly.graph_objects as go
import io
import re
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from Bio.Align import PairwiseAligner
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ============================================================
# ALLERGENLENS
# Interactive Bioinformatics Platform for Food-Allergen Analysis
# ============================================================

st.set_page_config(
    page_title="AllergenLens | Bioinformatics Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

# -------------------- SESSION STATE --------------------------

DEFAULTS = {
    "protein_records": [],
    "validated_sequence": "",
    "validated_name": "",
    "blast_rid": "",
    "blast_query_length": 0,
    "blast_query_name": "",
    "blast_submitted_at": None,
    "blast_last_checked": 0.0,
    "blast_status": "",
    "blast_results": [],
    "pairwise_results": None,
    "analysis_log": [],
}

for key, value in DEFAULTS.items():
    if key not in st.session_state:
        st.session_state[key] = value


def log_event(message):
    st.session_state.analysis_log.append({
        "Time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "Activity": message,
    })


# -------------------- GENERAL HELPERS -----------------------

def clean_sequence(sequence):
    """Remove FASTA headers, whitespace, and non-letter characters."""
    lines = [
        line.strip()
        for line in sequence.splitlines()
        if line.strip() and not line.startswith(">")
    ]
    raw = "".join(lines).upper()
    return re.sub(r"[^A-Z]", "", raw)


def validate_protein_sequence(sequence):
    """Validate a protein sequence using the standard amino-acid alphabet."""
    seq = clean_sequence(sequence)
    valid = set(AMINO_ACIDS)

    if not seq:
        return False, seq, "No sequence was provided."

    invalid = sorted(set(seq) - valid)

    if invalid:
        return (
            False,
            seq,
            "Unsupported character(s): " + ", ".join(invalid)
            + ". Use the 20 standard amino-acid letters only.",
        )

    return True, seq, "Sequence passed basic character validation."


def download_csv(df, filename):
    return df.to_csv(index=False).encode("utf-8")


def amino_acid_composition(sequence):
    seq = clean_sequence(sequence)
    length = len(seq)

    if length == 0:
        return pd.DataFrame(columns=["Amino Acid", "Count", "Percentage"])

    counts = {aa: seq.count(aa) for aa in AMINO_ACIDS}

    return pd.DataFrame({
        "Amino Acid": list(counts.keys()),
        "Count": list(counts.values()),
        "Percentage": [
            counts[aa] / length * 100 for aa in AMINO_ACIDS
        ],
    })


def xml_find_text(element, tag_name):
    """Find a descendant by local XML tag name, ignoring namespaces."""
    if element is None:
        return ""

    for node in element.iter():
        if node.tag.split("}")[-1].lower() == tag_name.lower():
            if node.text:
                return node.text.strip()

    return ""


# -------------------- UNIPROT HELPERS -----------------------

def search_uniprot(query, limit=10):
    params = {
        "query": query,
        "format": "json",
        "size": limit,
        "fields": (
            "accession,id,protein_name,organism_name,"
            "length,sequence,gene_names"
        ),
    }

    response = requests.get(
        UNIPROT_URL,
        params=params,
        timeout=30,
        headers={"User-Agent": "AllergenLens/1.0"},
    )
    response.raise_for_status()

    data = response.json()
    records = []

    for item in data.get("results", []):
        protein_description = (
            item.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", "")
        )

        if not protein_description:
            descriptions = (
                item.get("proteinDescription", {})
                .get("submissionNames", [])
            )
            if descriptions:
                protein_description = (
                    descriptions[0].get("fullName", {}).get("value", "")
                )

        organism = (
            item.get("organism", {}).get("scientificName", "Unknown")
        )
        sequence = item.get("sequence", {}).get("value", "")
        genes = item.get("genes", [])
        gene_name = ""

        if genes:
            gene_name = (
                genes[0].get("geneName", {}).get("value", "")
            )

        records.append({
            "Accession": item.get("primaryAccession", ""),
            "Entry Name": item.get("uniProtkbId", ""),
            "Protein": protein_description or "Protein name unavailable",
            "Organism": organism,
            "Gene": gene_name,
            "Length": len(sequence),
            "Sequence": sequence,
            "Source": "UniProt",
        })

    return records


# -------------------- PAIRWISE ALIGNMENT --------------------

def calculate_pairwise_alignment(sequence_a, sequence_b):
    """
    Run a global alignment and calculate alignment identity and coverage.
    Identity is calculated as identical aligned columns divided by total
    alignment columns, including gap columns.
    """
    seq_a = clean_sequence(sequence_a)
    seq_b = clean_sequence(sequence_b)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(seq_a, seq_b)[0]
    coordinates = alignment.coordinates

    aligned_a = []
    aligned_b = []

    for i in range(len(coordinates[0]) - 1):
        a_start = int(coordinates[0][i])
        a_end = int(coordinates[0][i + 1])
        b_start = int(coordinates[1][i])
        b_end = int(coordinates[1][i + 1])

        a_step = a_end - a_start
        b_step = b_end - b_start

        if a_step > 0 and b_step > 0:
            # Both sequences advance in an aligned segment.
            aligned_a.append(seq_a[a_start:a_end])
            aligned_b.append(seq_b[b_start:b_end])

        elif a_step > 0:
            aligned_a.append(seq_a[a_start:a_end])
            aligned_b.append("-" * a_step)

        elif b_step > 0:
            aligned_a.append("-" * b_step)
            aligned_b.append(seq_b[b_start:b_end])

    aligned_a = "".join(aligned_a)
    aligned_b = "".join(aligned_b)

    alignment_length = min(len(aligned_a), len(aligned_b))
    identical = sum(
        1 for a, b in zip(aligned_a, aligned_b)
        if a == b and a != "-" and b != "-"
    )

    identity = (
        identical / alignment_length * 100
        if alignment_length else 0
    )

    nongap_a = sum(1 for aa in aligned_a if aa != "-")
    nongap_b = sum(1 for aa in aligned_b if aa != "-")

    coverage_a = nongap_a / len(seq_a) * 100 if seq_a else 0
    coverage_b = nongap_b / len(seq_b) * 100 if seq_b else 0

    return {
        "score": float(alignment.score),
        "identity": identity,
        "alignment_length": alignment_length,
        "identical_residues": identical,
        "coverage_a": coverage_a,
        "coverage_b": coverage_b,
        "aligned_a": aligned_a,
        "aligned_b": aligned_b,
    }


# -------------------- NCBI BLASTp HELPERS -------------------

def submit_blast(sequence, database, email=""):
    params = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": database,
        "QUERY": ">AllergenLens_query\n" + clean_sequence(sequence),
        "FORMAT_TYPE": "XML2",
        "HITLIST_SIZE": 50,
        "EXPECT": 10,
        "TOOL": "AllergenLens",
    }

    if email.strip():
        params["EMAIL"] = email.strip()

    response = requests.post(
        BLAST_URL,
        data=params,
        timeout=60,
        headers={"User-Agent": "AllergenLens/1.0"},
    )
    response.raise_for_status()

    rid_match = re.search(r"RID\s*=\s*([A-Z0-9-]+)", response.text)
    rtoe_match = re.search(r"RTOE\s*=\s*(\d+)", response.text)

    if not rid_match:
        raise ValueError(
            "NCBI did not return a request ID. Please try again later "
            "or check the sequence and BLAST service."
        )

    return rid_match.group(1), int(rtoe_match.group(1)) if rtoe_match else None


def check_blast_status(rid):
    response = requests.get(
        BLAST_URL,
        params={
            "CMD": "Get",
            "RID": rid,
            "FORMAT_OBJECT": "SearchInfo",
        },
        timeout=30,
        headers={"User-Agent": "AllergenLens/1.0"},
    )
    response.raise_for_status()

    text = response.text.upper()

    if "STATUS=WAITING" in text or "STATUS = WAITING" in text:
        return "WAITING"
    if "STATUS=READY" in text or "STATUS = READY" in text:
        if "THEREARENOMESSAGES" in text:
            return "READY"
        return "READY"
    if "STATUS=FAILED" in text or "STATUS = FAILED" in text:
        return "FAILED"
    if "STATUS=UNKNOWN" in text or "STATUS = UNKNOWN" in text:
        return "UNKNOWN"

    # Handle responses where the status formatting differs.
    if "WAITING" in text:
        return "WAITING"
    if "FAILED" in text:
        return "FAILED"
    if "READY" in text:
        return "READY"

    return "UNKNOWN"


def fetch_blast_xml(rid):
    response = requests.get(
        BLAST_URL,
        params={
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "XML2",
        },
        timeout=60,
        headers={"User-Agent": "AllergenLens/1.0"},
    )
    response.raise_for_status()
    return response.text


def parse_blast_xml(xml_text):
    root = ET.fromstring(xml_text)
    rows = []

    for hit in root.iter():
        if hit.tag.split("}")[-1].lower() != "hit":
            continue

        accession = xml_find_text(hit, "accession")
        title = xml_find_text(hit, "title")

        if not title:
            title = xml_find_text(hit, "id")

        hsps = [
            node for node in hit.iter()
            if node.tag.split("}")[-1].lower() == "hsp"
        ]

        for hsp in hsps:
            identity = xml_find_text(hsp, "identity")
            align_len = xml_find_text(hsp, "align-len")
            query_from = xml_find_text(hsp, "query-from")
            query_to = xml_find_text(hsp, "query-to")
            evalue = xml_find_text(hsp, "evalue")
            bit_score = xml_find_text(hsp, "bit-score")

            try:
                ident_num = int(identity)
                align_num = int(align_len)
                query_from_num = int(query_from)
                query_to_num = int(query_to)
                evalue_num = float(evalue)
                bit_score_num = float(bit_score)

                identity_pct = (
                    ident_num / align_num * 100 if align_num else 0
                )
                query_span = abs(query_to_num - query_from_num) + 1
                query_coverage = (
                    query_span
                    / st.session_state.blast_query_length
                    * 100
                    if st.session_state.blast_query_length else 0
                )

            except (ValueError, TypeError, ZeroDivisionError):
                continue

            rows.append({
                "Hit": title or "Unknown hit",
                "Accession": accession,
                "Identity (%)": round(identity_pct, 2),
                "Query Coverage (%)": round(query_coverage, 2),
                "Alignment Length": align_num,
                "E-value": evalue_num,
                "Bit Score": round(bit_score_num, 2),
            })

    # Keep the best-scoring HSP for each accession/title.
    deduplicated = {}
    for row in rows:
        key = row["Accession"] or row["Hit"]
        if (
            key not in deduplicated
            or row["Bit Score"] > deduplicated[key]["Bit Score"]
        ):
            deduplicated[key] = row

    results = list(deduplicated.values())
    results.sort(key=lambda x: x["Bit Score"], reverse=True)
    return results


# -------------------- SIDEBAR -------------------------------

with st.sidebar:
    st.title("🧬 AllergenLens")
    st.caption("Food-allergen protein exploration")

    page = st.radio(
        "Navigation",
        [
            "Home",
            "Results Dashboard",
            "Database Explorer",
            "Sequence Laboratory",
            "BLASTp Analysis",
            "Protein Comparison",
            "Evidence & Limitations",
        ],
    )

    st.divider()
    st.caption("Bioinformatics learning and research prototype")
    st.caption("Results depend on input sequences and external services.")


# -------------------- HOME -----------------------------------

if page == "Home":
    st.title("🧬 AllergenLens")
    st.subheader(
        "An Interactive Bioinformatics Platform for Food-Allergen Analysis"
    )

    st.write(
        "Explore protein records, inspect amino-acid sequences, perform "
        "BLASTp similarity searches, compare two proteins, and visualize "
        "the results through one interactive interface."
    )

    c1, c2, c3 = st.columns(3)
    c1.metric("Core workflow", "Explore → Analyze")
    c2.metric("Sequence type", "Protein")
    c3.metric("Primary database", "UniProt")

    st.divider()
    st.markdown("### Explore the platform")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("#### 1. Database Explorer")
        st.write(
            "Search UniProt for protein records, inspect annotations, "
            "and download protein sequences."
        )

        st.markdown("#### 2. Sequence Laboratory")
        st.write(
            "Validate a protein sequence, examine amino-acid composition, "
            "and export a FASTA file."
        )

        st.markdown("#### 3. BLASTp Analysis")
        st.write(
            "Submit a protein sequence to NCBI BLASTp and examine returned "
            "hits, identity, query coverage, E-values, and bit scores."
        )

    with col2:
        st.markdown("#### 4. Protein Comparison")
        st.write(
            "Perform a global pairwise alignment and inspect identity, "
            "coverage, alignment score, and aligned sequences."
        )

        st.markdown("#### 5. Results Dashboard")
        st.write(
            "Review analyses saved in the current app session, compare "
            "protein lengths, and visualize BLASTp results."
        )

        st.markdown("#### 6. Evidence & Limitations")
        st.write(
            "Understand the difference between sequence similarity and "
            "clinical allergenicity or cross-reactivity."
        )

    st.info(
        "Important: sequence similarity alone cannot confirm that two "
        "proteins cause clinical cross-reactivity. Laboratory and clinical "
        "evidence are needed for such conclusions."
    )


# -------------------- RESULTS DASHBOARD ----------------------

elif page == "Results Dashboard":
    st.title("📊 Results Dashboard")
    st.write(
        "A live summary of the records and analyses saved in this "
        "current Streamlit session."
    )

    records = st.session_state.protein_records
    blast_df = pd.DataFrame(st.session_state.blast_results)
    pairwise = st.session_state.pairwise_results

    unique_accessions = {
        r.get("Accession", "")
        for r in records
        if r.get("Accession", "")
    }

    total_proteins = len(unique_accessions)
    total_blast_hits = len(blast_df) if not blast_df.empty else 0
    pairwise_count = 1 if pairwise else 0
    sequence_count = 1 if st.session_state.validated_sequence else 0

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Unique UniProt records", total_proteins)
    m2.metric("Validated sequence available", "Yes" if sequence_count else "No")
    m3.metric("Stored BLASTp hits", total_blast_hits)
    m4.metric("Pairwise comparisons", pairwise_count)

    st.divider()
    st.markdown("### Protein records")

    if records:
        protein_df = pd.DataFrame([
            {
                "Accession": r.get("Accession", ""),
                "Entry Name": r.get("Entry Name", ""),
                "Protein": r.get("Protein", ""),
                "Organism": r.get("Organism", ""),
                "Gene": r.get("Gene", ""),
                "Length": r.get("Length", 0),
            }
            for r in records
        ]).drop_duplicates(subset=["Accession"])

        st.dataframe(protein_df, use_container_width=True, hide_index=True)

        if not protein_df.empty:
            length_chart = px.bar(
                protein_df,
                x="Accession",
                y="Length",
                color="Organism",
                hover_data=["Protein", "Gene"],
                title="Protein sequence lengths",
                labels={"Length": "Length (amino acids)"},
            )
            st.plotly_chart(length_chart, use_container_width=True)

        st.download_button(
            "Download protein record table (CSV)",
            data=download_csv(protein_df, "allergenlens_protein_records.csv"),
            file_name="allergenlens_protein_records.csv",
            mime="text/csv",
        )
    else:
        st.info(
            "No UniProt records have been loaded in this session yet. "
            "Use Database Explorer to search for proteins."
        )

    st.divider()
    st.markdown("### BLASTp summary")

    if not blast_df.empty:
        st.dataframe(blast_df, use_container_width=True, hide_index=True)

        numeric_cols = [
            "Identity (%)",
            "Query Coverage (%)",
            "Bit Score",
        ]
        available_cols = [
            c for c in numeric_cols if c in blast_df.columns
        ]

        if available_cols:
            chart_col = st.selectbox(
                "Choose a BLASTp metric",
                available_cols,
                key="dashboard_blast_metric",
            )
            chart_df = blast_df.copy()
            chart_df["Short Hit"] = chart_df["Hit"].astype(str).str.slice(0, 45)

            fig = px.bar(
                chart_df.head(20),
                x="Short Hit",
                y=chart_col,
                hover_data=["Accession", "E-value"],
                title=f"Top BLASTp hits by {chart_col}",
            )
            fig.update_layout(xaxis_tickangle=-35)
            st.plotly_chart(fig, use_container_width=True)

        st.download_button(
            "Download BLASTp results (CSV)",
            data=download_csv(blast_df, "allergenlens_blastp_results.csv"),
            file_name="allergenlens_blastp_results.csv",
            mime="text/csv",
        )
    else:
        st.info(
            "No BLASTp results are stored in this session. Submit a search "
            "from the BLASTp Analysis page."
        )

    st.divider()
    st.markdown("### Latest pairwise comparison")

    if pairwise:
        a, b, c, d = st.columns(4)
        a.metric("Identity", f"{pairwise['identity']:.2f}%")
        b.metric("Alignment score", f"{pairwise['score']:.2f}")
        c.metric("Alignment length", pairwise["alignment_length"])
        d.metric("Identical residues", pairwise["identical_residues"])

        st.caption(
            "Identity uses identical aligned columns divided by the total "
            "alignment length, including gaps. The alignment score is not "
            "a percentage."
        )
    else:
        st.info("Run a comparison in Protein Comparison to populate this panel.")

    st.divider()
    st.markdown("### Session activity")

    if st.session_state.analysis_log:
        log_df = pd.DataFrame(st.session_state.analysis_log[::-1])
        st.dataframe(log_df, use_container_width=True, hide_index=True)

        st.download_button(
            "Download activity log (CSV)",
            data=download_csv(log_df, "allergenlens_activity_log.csv"),
            file_name="allergenlens_activity_log.csv",
            mime="text/csv",
        )
    else:
        st.caption("Activities will appear here as you use the platform.")


# -------------------- DATABASE EXPLORER ----------------------

elif page == "Database Explorer":
    st.title("🔎 UniProt Database Explorer")
    st.write(
        "Search the UniProtKB protein database. You can search by protein "
        "name, organism, gene, or accession."
    )

    with st.form("uniprot_search_form"):
        query = st.text_input(
            "Search query",
            value="",
            placeholder='Example: Ara h 1 OR organism:"Arachis hypogaea"',
        )
        limit = st.slider("Maximum records", 5, 50, 10, step=5)
        submitted = st.form_submit_button("Search UniProt")

    if submitted:
        if not query.strip():
            st.warning("Enter a search term first.")
        else:
            with st.spinner("Searching UniProt..."):
                try:
                    results = search_uniprot(query.strip(), limit)
                    st.session_state.protein_records = results

                    if results:
                        log_event(f"UniProt search: {query.strip()}")
                        st.success(f"Retrieved {len(results)} protein records.")
                    else:
                        st.warning(
                            "No records were returned. Try a broader search."
                        )

                except requests.RequestException as exc:
                    st.error(f"Could not reach UniProt: {exc}")
                except Exception as exc:
                    st.error(f"Search error: {exc}")

    records = st.session_state.protein_records

    if records:
        df = pd.DataFrame([
            {
                "Accession": r["Accession"],
                "Entry Name": r["Entry Name"],
                "Protein": r["Protein"],
                "Organism": r["Organism"],
                "Gene": r["Gene"],
                "Length": r["Length"],
            }
            for r in records
        ])

        st.markdown("### Retrieved records")
        st.dataframe(df, use_container_width=True, hide_index=True)

        selected_accession = st.selectbox(
            "Select a protein record",
            options=[r["Accession"] for r in records],
            format_func=lambda acc: next(
                (
                    f"{r['Accession']} — {r['Protein']} ({r['Organism']})"
                    for r in records if r["Accession"] == acc
                ),
                acc,
            ),
        )

        selected = next(
            r for r in records if r["Accession"] == selected_accession
        )

        st.markdown("### Selected protein")

        x1, x2, x3 = st.columns(3)
        x1.metric("Sequence length", selected["Length"])
        x2.metric("Organism", selected["Organism"])
        x3.metric("Accession", selected["Accession"])

        st.write("**Protein:**", selected["Protein"])
        st.write("**Entry name:**", selected["Entry Name"])
        st.write("**Gene:**", selected["Gene"] or "Not provided")

        st.text_area(
            "FASTA sequence",
            value=f">{selected['Accession']} | {selected['Protein']}\n"
                  f"{selected['Sequence']}",
            height=180,
            key="uniprot_selected_fasta",
        )

        st.download_button(
            "Download selected sequence (FASTA)",
            data=f">{selected['Accession']} | {selected['Protein']}\n"
                 f"{selected['Sequence']}\n",
            file_name=f"{selected['Accession']}.fasta",
            mime="text/plain",
        )

        if st.button("Use selected sequence in Sequence Laboratory"):
            st.session_state.validated_sequence = selected["Sequence"]
            st.session_state.validated_name = (
                f"{selected['Accession']} — {selected['Protein']}"
            )
            log_event(f"Selected UniProt sequence: {selected['Accession']}")
            st.success(
                "Sequence loaded. Open Sequence Laboratory or BLASTp Analysis."
            )

        st.download_button(
            "Download all displayed records (CSV)",
            data=download_csv(df, "allergenlens_uniprot_records.csv"),
            file_name="allergenlens_uniprot_records.csv",
            mime="text/csv",
        )


# -------------------- SEQUENCE LABORATORY --------------------

elif page == "Sequence Laboratory":
    st.title("🧪 Sequence Laboratory")
    st.write(
        "Validate protein sequence characters, inspect amino-acid "
        "composition, and export sequence data."
    )

    default_sequence = st.session_state.validated_sequence

    sequence_input = st.text_area(
        "Paste a protein sequence or FASTA record",
        value=default_sequence,
        height=200,
        placeholder=">Protein_name\nMKWVTFISLLFLFSSAYS...",
    )

    sequence_name = st.text_input(
        "Protein label",
        value=st.session_state.validated_name or "User-provided protein",
    )

    if st.button("Validate and analyze sequence", type="primary"):
        valid, seq, message = validate_protein_sequence(sequence_input)

        if valid:
            st.session_state.validated_sequence = seq
            st.session_state.validated_name = sequence_name.strip()
            log_event(f"Validated sequence: {sequence_name.strip()}")
            st.success(message)
        else:
            st.error(message)

    seq = st.session_state.validated_sequence

    if seq:
        st.markdown("### Sequence summary")

        a, b, c = st.columns(3)
        a.metric("Sequence length", len(seq))
        b.metric("Distinct amino acids", len(set(seq)))
        c.metric(
            "Most frequent amino acid",
            max(AMINO_ACIDS, key=lambda aa: seq.count(aa)),
        )

        comp_df = amino_acid_composition(seq)

        fig = px.bar(
            comp_df,
            x="Amino Acid",
            y="Percentage",
            hover_data=["Count"],
            title="Amino-acid composition",
            labels={"Percentage": "Composition (%)"},
        )
        st.plotly_chart(fig, use_container_width=True)

        st.dataframe(comp_df, use_container_width=True, hide_index=True)

        st.download_button(
            "Download composition table (CSV)",
            data=download_csv(comp_df, "amino_acid_composition.csv"),
            file_name="amino_acid_composition.csv",
            mime="text/csv",
        )

        st.download_button(
            "Download validated sequence (FASTA)",
            data=f">{st.session_state.validated_name or 'protein'}\n{seq}\n",
            file_name="validated_protein.fasta",
            mime="text/plain",
        )

        with st.expander("View normalized sequence"):
            st.code(seq, language=None)

        st.caption(
            "This page performs basic sequence-character validation and "
            "composition analysis. It does not determine allergenicity."
        )


# -------------------- BLASTp ANALYSIS ------------------------

elif page == "BLASTp Analysis":
    st.title("🧬 BLASTp Analysis")
    st.write(
        "Compare a protein sequence against a protein database using "
        "the NCBI BLAST service."
    )

    st.warning(
        "NCBI BLAST is an external service. Searches may take time, and "
        "availability or result formatting can vary. Avoid submitting "
        "sensitive or unpublished sequences unless you are comfortable "
        "with the service's data policies."
    )

    if st.session_state.validated_sequence:
        st.caption(
            "A validated sequence is available: "
            f"{st.session_state.validated_name or 'Unnamed protein'} "
            f"({len(st.session_state.validated_sequence)} amino acids)."
        )

    input_mode = st.radio(
        "Sequence source",
        ["Use validated sequence", "Paste a sequence"],
        horizontal=True,
    )

    if input_mode == "Use validated sequence":
        blast_sequence = st.session_state.validated_sequence
        if blast_sequence:
            st.text_area(
                "Sequence preview",
                value=blast_sequence[:1000],
                height=120,
                disabled=True,
            )
        else:
            st.info(
                "No validated sequence is available. Load one from "
                "Database Explorer or Sequence Laboratory."
            )
    else:
        blast_sequence = st.text_area(
            "Protein sequence",
            height=160,
            placeholder=">Query\nMKWVTFISLLFLFSSAYS...",
            key="blast_manual_sequence",
        )

    db_choice = st.selectbox(
        "BLAST database",
        options=["nr", "swissprot"],
        format_func=lambda x: (
            "NCBI nr — broad protein database"
            if x == "nr"
            else "Swiss-Prot — curated protein database"
        ),
    )

    email = st.text_input(
        "Email (optional)",
        placeholder="you@example.com",
        help="An email may help identify the request. Follow NCBI's usage guidance.",
    )

    if st.button("Submit BLASTp search", type="primary"):
        valid, seq, message = validate_protein_sequence(blast_sequence or "")

        if not valid:
            st.error(message)
        elif len(seq) < 5:
            st.error("The sequence is very short. Check the input.")
        else:
            try:
                with st.spinner("Submitting request to NCBI..."):
                    rid, rtoe = submit_blast(seq, db_choice, email)

                st.session_state.blast_rid = rid
                st.session_state.blast_query_length = len(seq)
                st.session_state.blast_query_name = (
                    st.session_state.validated_name
                    if input_mode == "Use validated sequence"
                    else "Manually entered sequence"
                )
                st.session_state.blast_submitted_at = datetime.now(
                    timezone.utc
                ).isoformat()
                st.session_state.blast_last_checked = 0.0
                st.session_state.blast_status = "WAITING"
                st.session_state.blast_results = []

                log_event(f"Submitted BLASTp search ({db_choice}), RID {rid}")
                st.success(f"BLASTp request submitted. Request ID: {rid}")

                if rtoe:
                    st.caption(
                        f"NCBI's estimated processing time: about {rtoe} seconds."
                    )

            except Exception as exc:
                st.error(f"BLASTp submission failed: {exc}")

    if st.session_state.blast_rid:
        st.divider()
        st.markdown("### Current BLASTp request")

        st.write("**Request ID:**", st.session_state.blast_rid)
        st.write("**Query:**", st.session_state.blast_query_name)
        st.write(
            "**Submitted:**",
            st.session_state.blast_submitted_at or "Unknown",
        )
        st.write("**Status:**", st.session_state.blast_status or "Not checked")

        elapsed = (
            time.time() - st.session_state.blast_last_checked
            if st.session_state.blast_last_checked else None
        )

        if elapsed is not None and elapsed < 60:
            st.info(
                f"Please wait approximately {int(60 - elapsed)} more seconds "
                "before checking the request again, to avoid excessive polling."
            )

        if st.button("Check BLASTp status / retrieve results"):
            elapsed = (
                time.time() - st.session_state.blast_last_checked
                if st.session_state.blast_last_checked else 999
            )

            if elapsed < 60:
                st.warning(
                    "Please wait at least 60 seconds between status checks."
                )
            else:
                try:
                    with st.spinner("Checking NCBI request..."):
                        status = check_blast_status(
                            st.session_state.blast_rid
                        )

                    st.session_state.blast_last_checked = time.time()
                    st.session_state.blast_status = status

                    if status == "WAITING":
                        st.info(
                            "NCBI is still processing the request. Wait and "
                            "check again later."
                        )

                    elif status == "READY":
                        with st.spinner("Retrieving BLASTp results..."):
                            xml_result = fetch_blast_xml(
                                st.session_state.blast_rid
                            )
                            hits = parse_blast_xml(xml_result)

                        st.session_state.blast_results = hits
                        log_event(
                            f"Retrieved BLASTp results: {len(hits)} hits"
                        )

                        if hits:
                            st.success(f"Retrieved {len(hits)} BLASTp hits.")
                        else:
                            st.warning(
                                "The request is ready, but no hits were parsed. "
                                "The returned XML format may have changed, or "
                                "there may be no reportable matches. Try again "
                                "or inspect the request on the NCBI BLAST website."
                            )

                    elif status == "FAILED":
                        st.error(
                            "NCBI reported that the BLAST request failed. "
                            "Try submitting it again later."
                        )

                    else:
                        st.warning(
                            "The request status could not be interpreted. "
                            "Wait a little and retry."
                        )

                except Exception as exc:
                    st.error(f"Could not check/retrieve BLASTp results: {exc}")

        blast_results = st.session_state.blast_results

        if blast_results:
            blast_df = pd.DataFrame(blast_results)

            st.markdown("### BLASTp hits")
            st.dataframe(blast_df, use_container_width=True, hide_index=True)

            fig = px.scatter(
                blast_df,
                x="Query Coverage (%)",
                y="Identity (%)",
                size="Bit Score",
                hover_name="Hit",
                hover_data=["Accession", "E-value", "Alignment Length"],
                title="BLASTp hit identity versus query coverage",
            )
            st.plotly_chart(fig, use_container_width=True)

            st.download_button(
                "Download BLASTp results (CSV)",
                data=download_csv(blast_df, "blastp_results.csv"),
                file_name="blastp_results.csv",
                mime="text/csv",
            )

            st.caption(
                "Identity is the percentage of identical residues in an "
                "aligned region. Query coverage estimates how much of the "
                "query is covered by the reported HSP. E-values estimate "
                "the number of alignments of similar quality expected by "
                "chance in the searched database. Interpret these metrics "
                "together; they do not establish clinical cross-reactivity."
            )


# -------------------- PROTEIN COMPARISON ---------------------

elif page == "Protein Comparison":
    st.title("🔬 Protein Comparison")
    st.write(
        "Compare two protein sequences using global pairwise alignment."
    )

    seq_options = ["Enter manually", "Use validated sequence"]

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("### Protein A")
        source_a = st.selectbox(
            "Sequence source A",
            seq_options,
            key="source_a",
        )

        if source_a == "Use validated sequence":
            seq_a_input = st.text_area(
                "Protein A sequence",
                value=st.session_state.validated_sequence,
                height=160,
                key="seq_a_validated",
            )
        else:
            seq_a_input = st.text_area(
                "Protein A sequence",
                height=160,
                key="seq_a_manual",
                placeholder="Paste sequence A",
            )

        name_a = st.text_input("Protein A label", value="Protein A")

    with col2:
        st.markdown("### Protein B")
        seq_b_input = st.text_area(
            "Protein B sequence",
            height=160,
            key="seq_b_manual",
            placeholder="Paste sequence B",
        )
        name_b = st.text_input("Protein B label", value="Protein B")

    if st.button("Compare proteins", type="primary"):
        valid_a, seq_a, message_a = validate_protein_sequence(seq_a_input or "")
        valid_b, seq_b, message_b = validate_protein_sequence(seq_b_input or "")

        if not valid_a:
            st.error("Protein A: " + message_a)
        elif not valid_b:
            st.error("Protein B: " + message_b)
        else:
            with st.spinner("Aligning sequences..."):
                result = calculate_pairwise_alignment(seq_a, seq_b)

            result["name_a"] = name_a
            result["name_b"] = name_b
            result["length_a"] = len(seq_a)
            result["length_b"] = len(seq_b)
            st.session_state.pairwise_results = result
            log_event(f"Pairwise comparison: {name_a} vs {name_b}")
            st.success("Pairwise alignment completed.")

    result = st.session_state.pairwise_results

    if result:
        st.divider()
        st.markdown("### Comparison summary")

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Identity", f"{result['identity']:.2f}%")
        c2.metric("Alignment score", f"{result['score']:.2f}")
        c3.metric("Alignment length", result["alignment_length"])
        c4.metric("Identical residues", result["identical_residues"])

        coverage_df = pd.DataFrame({
            "Protein": [result["name_a"], result["name_b"]],
            "Sequence length": [result["length_a"], result["length_b"]],
            "Coverage (%)": [
                result["coverage_a"], result["coverage_b"]
            ],
        })

        st.markdown("### Sequence coverage")
        fig = px.bar(
            coverage_df,
            x="Protein",
            y="Coverage (%)",
            hover_data=["Sequence length"],
            title="Coverage of each protein in the global alignment",
            range_y=[0, 100],
        )
        st.plotly_chart(fig, use_container_width=True)

        st.dataframe(coverage_df, use_container_width=True, hide_index=True)

        st.markdown("### Aligned sequences")
        st.caption(
            "Dashes represent alignment gaps. Long sequences are displayed "
            "in wrapped blocks for readability."
        )

        chunk_size = 70
        aligned_a = result["aligned_a"]
        aligned_b = result["aligned_b"]

        for start in range(0, len(aligned_a), chunk_size):
            end = start + chunk_size
            st.code(
                f"{result['name_a']}: {aligned_a[start:end]}\n"
                f"{result['name_b']}: {aligned_b[start:end]}",
                language=None,
            )

        comparison_export = pd.DataFrame([{
            "Protein A": result["name_a"],
            "Protein B": result["name_b"],
            "Length A": result["length_a"],
            "Length B": result["length_b"],
            "Identity (%)": round(result["identity"], 2),
            "Alignment score": round(result["score"], 2),
            "Alignment length": result["alignment_length"],
            "Identical residues": result["identical_residues"],
            "Coverage A (%)": round(result["coverage_a"], 2),
            "Coverage B (%)": round(result["coverage_b"], 2),
        }])

        st.download_button(
            "Download comparison summary (CSV)",
            data=download_csv(comparison_export, "pairwise_comparison.csv"),
            file_name="pairwise_comparison.csv",
            mime="text/csv",
        )

        st.caption(
            "This is a global pairwise alignment using a simple scoring "
            "scheme. Its score depends on the scoring parameters and is "
            "not directly comparable with a BLAST bit score. Sequence "
            "similarity alone cannot confirm allergenicity or clinical "
            "cross-reactivity."
        )


# -------------------- EVIDENCE & LIMITATIONS -----------------

elif page == "Evidence & Limitations":
    st.title("📚 Evidence & Limitations")
    st.write(
        "This section explains what AllergenLens results can and cannot "
        "establish."
    )

    st.markdown("### Data sources and methods")

    st.markdown("""
    - **UniProtKB:** protein records, accessions, annotations, and sequences.
    - **NCBI BLASTp:** similarity searching against a selected protein database.
    - **Biopython PairwiseAligner:** global pairwise sequence alignment.
    - **Pandas:** tabular result processing and CSV exports.
    - **Plotly:** interactive data visualizations.
    - **Streamlit:** interactive application interface.
    """)

    st.markdown("### Interpreting common metrics")

    metric_df = pd.DataFrame([
        {
            "Metric": "Sequence identity",
            "Meaning": "Percentage of identical residues within an alignment.",
            "Caution": "Depends on the aligned region and alignment definition.",
        },
        {
            "Metric": "Query coverage",
            "Meaning": "Approximate proportion of the query covered by an alignment.",
            "Caution": "High coverage alone does not prove biological function.",
        },
        {
            "Metric": "E-value",
            "Meaning": "Expected number of chance matches of similar quality in the database search.",
            "Caution": "Depends on database size, scoring, and search settings.",
        },
        {
            "Metric": "Bit score",
            "Meaning": "Normalized alignment score reported by BLAST.",
            "Caution": "Do not treat it as a percentage identity.",
        },
        {
            "Metric": "Global alignment score",
            "Meaning": "Score produced by the pairwise alignment algorithm.",
            "Caution": "Depends on match, mismatch, and gap penalties.",
        },
    ])

    st.dataframe(metric_df, use_container_width=True, hide_index=True)

    st.markdown("### Important scientific limitations")

    st.markdown("""
    1. Protein sequence similarity is not equivalent to confirmed clinical
       allergenicity or cross-reactivity.
    2. Results depend on sequence quality, database selection, and tool settings.
    3. A negative BLASTp result does not prove that proteins have no shared
       structural or immunological properties.
    4. This platform is an exploratory educational tool, not a clinical
       diagnostic system or a substitute for allergy testing.
    5. A meaningful cross-reactivity assessment requires additional evidence,
       such as validated allergen annotations, epitope information, structural
       context, and relevant experimental or clinical studies.
    6. Session-based metrics and results may be lost when the app session ends.
    """)

    st.markdown("### Suggested references")

    st.markdown("""
    - UniProt: https://www.uniprot.org/
    - UniProt REST API documentation: https://www.uniprot.org/help/api
    - NCBI BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi
    - NCBI BLAST URL API guide: https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html
    - Biopython documentation: https://biopython.org/
    """)

    st.info(
        "For your project report, document the access date, accession numbers, "
        "database used, search parameters, and the actual outputs obtained. "
        "Do not report example or missing results as experimental findings."
    )


# -------------------- FOOTER ---------------------------------

st.divider()
st.caption(
    "AllergenLens | Educational bioinformatics prototype | "
    "External database and BLAST results depend on service availability."
)
