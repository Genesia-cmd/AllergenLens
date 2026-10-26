
import streamlit as st
import pandas as pd
import requests
import re
import time
import xml.etree.ElementTree as ET
from collections import Counter

st.set_page_config(
    page_title="AllergenLens",
    page_icon="🧬",
    layout="wide"
)

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

st.title("🧬 AllergenLens")
st.caption(
    "Interactive Food-Allergen Protein Exploration "
    "and Sequence Analysis Platform"
)

page = st.sidebar.radio(
    "Navigation",
    [
        "Home",
        "Database Explorer",
        "Sequence Laboratory",
        "BLASTp Analysis",
        "Protein Comparison",
        "Results & Visualization",
        "Evidence & Limitations"
    ]
)


# ==================================================
# HELPER FUNCTIONS
# ==================================================

@st.cache_data(ttl=3600)
def search_uniprot(query, limit):
    response = requests.get(
        UNIPROT_URL,
        params={
            "query": query,
            "format": "json",
            "size": limit,
            "fields": (
                "accession,protein_name,organism_name,"
                "gene_names,length,sequence"
            )
        },
        timeout=45
    )
    response.raise_for_status()
    return response.json().get("results", [])


def get_protein_name(record):
    description = record.get("proteinDescription", {})
    recommended = description.get("recommendedName", {})

    if recommended.get("fullName", {}).get("value"):
        return recommended["fullName"]["value"]

    alternatives = description.get("submissionNames", [])

    if alternatives:
        return alternatives[0].get("fullName", {}).get(
            "value", "Name unavailable"
        )

    return "Name unavailable"


def validate_protein(sequence_text):
    lines = [
        line.strip()
        for line in sequence_text.splitlines()
        if line.strip()
    ]

    headers = [
        line for line in lines
        if line.startswith(">")
    ]

    if len(headers) > 1:
        return None, "Please provide only one protein sequence."

    sequence = "".join(
        line for line in lines
        if not line.startswith(">")
    )

    sequence = re.sub(r"\s+", "", sequence).upper()

    if sequence.endswith("*"):
        sequence = sequence[:-1]

    allowed = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")

    if not sequence:
        return None, "No protein sequence was provided."

    invalid = sorted(set(sequence) - allowed)

    if invalid:
        return None, (
            "Invalid protein symbols: " + ", ".join(invalid)
        )

    if len(sequence) < 10:
        return None, (
            "Sequence is too short for this application. "
            "Please provide at least 10 amino acids."
        )

    return sequence, None


def xml_text(element, tag):
    """Find a child by tag name, ignoring XML namespaces."""
    for child in element.iter():
        if child.tag.split("}")[-1] == tag:
            return child.text or ""
    return ""


def parse_blast_xml(xml_content, query_length):
    root = ET.fromstring(xml_content)
    rows = []

    for hit in root.iter():
        if hit.tag.split("}")[-1] != "hit":
            continue

        descriptions = [
            element for element in hit.iter()
            if element.tag.split("}")[-1] == "description"
        ]

        title = ""
        accession = ""

        if descriptions:
            title = xml_text(descriptions[0], "title")
            accession = xml_text(descriptions[0], "accession")

        if not accession:
            accession = xml_text(hit, "accession")

        hsps = [
            element for element in hit.iter()
            if element.tag.split("}")[-1] == "hsp"
        ]

        if not hsps:
            continue

        # Report the best-scoring HSP for each hit.
        def bit_score(hsp):
            try:
                return float(xml_text(hsp, "bit-score"))
            except (ValueError, TypeError):
                return 0.0

        best = max(hsps, key=bit_score)

        try:
            identity = int(xml_text(best, "identity"))
            alignment_length = int(
                xml_text(best, "align-len")
            )
            query_start = int(
                xml_text(best, "query-from")
            )
            query_end = int(
                xml_text(best, "query-to")
            )
            evalue = float(xml_text(best, "evalue"))
            score = float(xml_text(best, "bit-score"))
        except (ValueError, TypeError):
            continue

        percent_identity = (
            100 * identity / alignment_length
            if alignment_length else 0
        )

        coverage = (
            100 * (abs(query_end - query_start) + 1)
            / query_length
            if query_length else 0
        )

        rows.append({
            "Protein / Hit": title or accession or "Unknown",
            "Accession": accession or "Unavailable",
            "Identity (%)": round(percent_identity, 2),
            "Query Coverage (%)": round(
                min(coverage, 100), 2
            ),
            "Alignment Length": alignment_length,
            "E-value": evalue,
            "Bit Score": round(score, 2)
        })

    return rows


def submit_blast(sequence, database, email):
    fasta = f">AllergenLens_query\n{sequence}"

    params = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": database,
        "QUERY": fasta,
        "TOOL": "AllergenLens",
        "FORMAT_TYPE": "XML2",
        "HITLIST_SIZE": 50,
        "EXPECT": 10
    }

    if email.strip():
        params["EMAIL"] = email.strip()

    response = requests.post(
        BLAST_URL,
        data=params,
        timeout=60
    )
    response.raise_for_status()

    rid_match = re.search(
        r"RID\s*=\s*([A-Z0-9-]+)",
        response.text
    )

    if not rid_match:
        raise ValueError(
            "NCBI did not return a request ID. "
            "Try again later or use the NCBI BLAST website."
        )

    rtoe_match = re.search(
        r"RTOE\s*=\s*(\d+)",
        response.text
    )

    return (
        rid_match.group(1),
        int(rtoe_match.group(1)) if rtoe_match else None
    )


def check_blast_status(rid):
    response = requests.get(
        BLAST_URL,
        params={
            "CMD": "Get",
            "RID": rid,
            "FORMAT_OBJECT": "SearchInfo"
        },
        timeout=45
    )
    response.raise_for_status()

    status_match = re.search(
        r"Status\s*=\s*(READY|WAITING|FAILED|UNKNOWN)",
        response.text,
        flags=re.IGNORECASE
    )

    if status_match:
        return status_match.group(1).upper()

    if "ThereAreHits=yes" in response.text:
        return "READY"

    return "UNKNOWN"


def retrieve_blast_results(rid):
    response = requests.get(
        BLAST_URL,
        params={
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "XML2"
        },
        timeout=90
    )
    response.raise_for_status()

    content = response.text

    if not content.strip():
        raise ValueError("NCBI returned an empty result.")

    return content


# ==================================================
# HOME
# ==================================================

if page == "Home":
    st.header("Welcome to AllergenLens")

    st.write(
        "Explore food-allergen proteins, retrieve biological "
        "sequences, perform similarity searches, compare "
        "proteins, and visualize computational results."
    )

    c1, c2, c3 = st.columns(3)
    c1.metric("Core Modules", "7")
    c2.metric("Protein Database", "UniProt")
    c3.metric("Similarity Search", "NCBI BLASTp")

    st.subheader("Platform Features")

    st.markdown("""
    - Search and explore protein records
    - Retrieve FASTA sequences
    - Validate protein sequence input
    - Run BLASTp similarity searches
    - Compare protein sequences
    - Visualize and export results
    """)


# ==================================================
# DATABASE EXPLORER
# ==================================================

elif page == "Database Explorer":
    st.header("🔎 Allergen Database Explorer")

    st.write(
        "Search UniProtKB by protein name, organism, "
        "or accession."
    )

    query = st.text_input(
        "Search proteins",
        placeholder="e.g. Ara h 2, peanut, P02754"
    )

    limit = st.selectbox(
        "Maximum records",
        [5, 10, 20, 50],
        index=1
    )

    if st.button("Search UniProt", type="primary"):
        if not query.strip():
            st.warning("Enter a search term.")
        else:
            try:
                with st.spinner("Searching UniProtKB..."):
                    records = search_uniprot(
                        query.strip(), limit
                    )

                rows = []

                for record in records:
                    seq_data = record.get("sequence", {})

                    rows.append({
                        "Accession": record.get(
                            "primaryAccession", ""
                        ),
                        "Protein": get_protein_name(record),
                        "Organism": record.get(
                            "organism", {}
                        ).get("scientificName", "Unknown"),
                        "Length (aa)": seq_data.get("length"),
                        "Reviewed": record.get(
                            "entryType", ""
                        ) == "UniProtKB reviewed (Swiss-Prot)"
                    })

                st.session_state["uniprot_records"] = records
                st.session_state["uniprot_table"] = pd.DataFrame(rows)

                if rows:
                    st.success(f"Retrieved {len(rows)} records.")
                else:
                    st.info("No records found.")

            except requests.RequestException as error:
                st.error(f"UniProt connection error: {error}")

    if "uniprot_table" in st.session_state:
        table = st.session_state["uniprot_table"]
        records = st.session_state["uniprot_records"]

        st.dataframe(
            table,
            use_container_width=True,
            hide_index=True
        )

        st.download_button(
            "Download results CSV",
            table.to_csv(index=False),
            "allergenlens_proteins.csv",
            "text/csv"
        )

        if records:
            options = {
                f"{r.get('primaryAccession', '')} — "
                f"{get_protein_name(r)}": i
                for i, r in enumerate(records)
            }

            selected = st.selectbox(
                "Choose a protein record",
                list(options.keys())
            )

            record = records[options[selected]]
            accession = record.get("primaryAccession", "")
            sequence = record.get("sequence", {}).get(
                "value", ""
            )

            st.write("**Organism:**", record.get(
                "organism", {}
            ).get("scientificName", "Unknown"))

            st.write("**Accession:**", accession)
            st.write("**Sequence length:**", len(sequence), "aa")

            if sequence:
                fasta = (
                    f">{accession} {get_protein_name(record)}\n"
                    + "\n".join(
                        sequence[i:i + 60]
                        for i in range(0, len(sequence), 60)
                    )
                )

                with st.expander("View FASTA"):
                    st.code(fasta, language=None)

                st.download_button(
                    "Download FASTA",
                    fasta,
                    file_name=f"{accession}.fasta",
                    mime="text/plain",
                    key="uniprot_fasta_download"
                )

                if st.button("Use this sequence in BLASTp"):
                    st.session_state["blast_input"] = fasta
                    st.success(
                        "Sequence prepared. Open BLASTp Analysis "
                        "in the sidebar."
                    )

            st.markdown(
                f"[Open UniProt record]"
                f"(https://www.uniprot.org/uniprotkb/{accession}/entry)"
            )


# ==================================================
# SEQUENCE LABORATORY
# ==================================================

elif page == "Sequence Laboratory":
    st.header("🧪 Sequence Analysis Laboratory")

    uploaded = st.file_uploader(
        "Upload a protein FASTA file",
        type=["fasta", "fa", "faa", "txt"]
    )

    pasted = st.text_area(
        "Or paste a protein sequence / FASTA",
        value=st.session_state.get("blast_input", ""),
        height=160,
        placeholder=">Protein_name\nMKWVTFISLL..."
    )

    raw = pasted

    if uploaded is not None:
        raw = uploaded.getvalue().decode(
            "utf-8", errors="replace"
        )

    if st.button("Validate Sequence", type="primary"):
        sequence, error = validate_protein(raw)

        if error:
            st.error(error)
        else:
            st.session_state["validated_sequence"] = sequence

            st.success("Sequence validation passed.")

            c1, c2 = st.columns(2)
            c1.metric("Length", f"{len(sequence)} aa")
            c2.metric("Unique residues", len(set(sequence)))

            counts = Counter(sequence)

            chart = pd.DataFrame({
                "Residue": list(counts.keys()),
                "Count": list(counts.values())
            }).sort_values("Residue")

            st.subheader("Amino-Acid Composition")
            st.bar_chart(chart.set_index("Residue")["Count"])

            fasta = ">AllergenLens_query\n" + "\n".join(
                sequence[i:i + 60]
                for i in range(0, len(sequence), 60)
            )

            st.download_button(
                "Download validated FASTA",
                fasta,
                file_name="allergenlens_query.fasta",
                mime="text/plain"
            )

            st.info(
                "Validation checks sequence formatting and symbols; "
                "it does not confirm allergen status."
            )


# ==================================================
# REAL NCBI BLASTp ANALYSIS
# ==================================================

elif page == "BLASTp Analysis":
    st.header("🧬 NCBI BLASTp Analysis")

    st.write(
        "Submit a protein sequence to NCBI and retrieve "
        "actual protein similarity-search results."
    )

    st.caption(
        "NCBI is a shared public service. Searches may take "
        "several minutes. Avoid repeatedly submitting the same query."
    )

    sequence_source = st.radio(
        "Sequence input",
        [
            "Use validated sequence",
            "Paste or upload a sequence"
        ]
    )

    sequence_text = ""

    if sequence_source == "Use validated sequence":
        if "validated_sequence" in st.session_state:
            sequence_text = st.session_state["validated_sequence"]
            st.success(
                f"Validated sequence loaded: {len(sequence_text)} aa"
            )
        else:
            st.info(
                "First validate a sequence in Sequence Laboratory, "
                "or select the paste/upload option."
            )

    if sequence_source == "Paste or upload a sequence":
        blast_upload = st.file_uploader(
            "Upload FASTA for BLASTp",
            type=["fasta", "fa", "faa", "txt"],
            key="blast_upload"
        )

        blast_paste = st.text_area(
            "Paste FASTA or amino-acid sequence",
            height=140,
            key="blast_paste"
        )

        sequence_text = blast_paste

        if blast_upload is not None:
            sequence_text = blast_upload.getvalue().decode(
                "utf-8", errors="replace"
            )

    database = st.selectbox(
        "NCBI protein database",
        ["nr", "swissprot"],
        format_func=lambda x: {
            "nr": "nr — broad protein database",
            "swissprot": "Swiss-Prot — reviewed protein records"
        }[x]
    )

    email = st.text_input(
        "Contact email (recommended by NCBI)",
        placeholder="you@example.com"
    )

    sequence, validation_error = validate_protein(sequence_text)

    if sequence is not None:
        st.caption(
            f"Input contains {len(sequence)} amino acids."
        )

    if validation_error and sequence_text.strip():
        st.warning(validation_error)

    if st.button("Submit BLASTp Search", type="primary"):
        if validation_error:
            st.error(validation_error)
        else:
            try:
                with st.spinner(
                    "Submitting your sequence to NCBI..."
                ):
                    rid, estimated_seconds = submit_blast(
                        sequence, database, email
                    )

                st.session_state["blast_rid"] = rid
                st.session_state["blast_query_length"] = len(sequence)
                st.session_state["blast_database"] = database
                st.session_state["blast_last_check"] = 0.0
                st.session_state.pop("blast_results", None)

                st.success("BLASTp job submitted successfully.")
                st.code(rid, language=None)

                st.write("**Database:**", database)
                st.write("**Query length:**", len(sequence), "aa")

                if estimated_seconds is not None:
                    st.write(
                        "NCBI's estimated processing time:",
                        f"{estimated_seconds} seconds"
                    )

                st.info(
                    "Save the request ID. Once processing has had "
                    "time to complete, use Check BLASTp Status. "
                    "NCBI recommends polling no more than once "
                    "per minute for an individual request."
                )

            except (requests.RequestException, ValueError) as error:
                st.error(f"BLASTp submission failed: {error}")

    if "blast_rid" in st.session_state:
        rid = st.session_state["blast_rid"]

        st.subheader("BLASTp Job")
        st.write("**Request ID:**", rid)
        st.write(
            "**Database:**",
            st.session_state.get("blast_database", "nr")
        )

        if st.button("Check BLASTp Status"):
            now = time.time()
            last_check = st.session_state.get(
                "blast_last_check", 0.0
            )

            if now - last_check < 60:
                remaining = int(60 - (now - last_check))
                st.warning(
                    f"Please wait approximately {remaining} seconds "
                    "before checking this job again."
                )
            else:
                st.session_state["blast_last_check"] = now

                try:
                    with st.spinner("Checking NCBI job status..."):
                        status = check_blast_status(rid)

                    st.write("**NCBI status:**", status)

                    if status == "WAITING":
                        st.info(
                            "Your search is still processing. "
                            "Wait at least one minute before checking again."
                        )

                    elif status == "FAILED":
                        st.error(
                            "NCBI reported that this search failed. "
                            "Try a new search later."
                        )

                    elif status == "UNKNOWN":
                        st.warning(
                            "NCBI returned an unknown status. "
                            "The job may not be available yet."
                        )

                    elif status == "READY":
                        with st.spinner("Retrieving BLASTp results..."):
                            xml_content = retrieve_blast_results(rid)

                        query_length = st.session_state[
                            "blast_query_length"
                        ]

                        rows = parse_blast_xml(
                            xml_content, query_length
                        )

                        st.session_state["blast_results"] = rows

                        if rows:
                            st.success(
                                f"Retrieved {len(rows)} protein hits."
                            )
                        else:
                            st.info(
                                "No parsed hits were found. "
                                "The query may have no reported matches "
                                "or the result format may have changed."
                            )

                except (
                    requests.RequestException,
                    ValueError,
                    ET.ParseError
                ) as error:
                    st.error(
                        f"Could not retrieve BLASTp results: {error}"
                    )

        if "blast_results" in st.session_state:
            results = st.session_state["blast_results"]

            if results:
                df = pd.DataFrame(results)

                st.subheader("Similarity Search Results")

                st.dataframe(
                    df,
                    use_container_width=True,
                    hide_index=True
                )

                st.download_button(
                    "Download BLASTp results CSV",
                    df.to_csv(index=False),
                    file_name="allergenlens_blastp_results.csv",
                    mime="text/csv"
                )

                st.subheader("Top Hits by Bit Score")

                chart_df = df.head(10).copy()
                chart_df["Display name"] = (
                    chart_df["Accession"].astype(str)
                )

                st.bar_chart(
                    chart_df.set_index("Display name")["Bit Score"]
                )

                st.caption(
                    "Identity is calculated from the best-scoring "
                    "high-scoring segment pair (HSP) per hit. "
                    "Query coverage here is the span of that HSP, "
                    "not necessarily the union of all aligned regions."
                )

                st.markdown(
                    f"[Open the NCBI BLAST results page]"
                    f"(https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID={rid})"
                )

            else:
                st.info("No hit table is currently available.")


# ==================================================
# PROTEIN COMPARISON
# ==================================================

elif page == "Protein Comparison":
    st.header("🧬 Protein Comparison")

    st.write(
        "Compare two protein sequences using local pairwise "
        "alignment. This is distinct from searching the NCBI database."
    )

    seq1_text = st.text_area(
        "First protein sequence",
        key="pair_seq1",
        height=100
    )

    seq2_text = st.text_area(
        "Second protein sequence",
        key="pair_seq2",
        height=100
    )

    if st.button("Align Protein Sequences"):
        seq1, error1 = validate_protein(seq1_text)
        seq2, error2 = validate_protein(seq2_text)

        if error1:
            st.error("First sequence: " + error1)
        elif error2:
            st.error("Second sequence: " + error2)
        else:
            try:
                from Bio.Align import PairwiseAligner

                aligner = PairwiseAligner()
                aligner.mode = "global"
                aligner.match_score = 2
                aligner.mismatch_score = -1
                aligner.open_gap_score = -5
                aligner.extend_gap_score = -0.5

                alignment = aligner.align(seq1, seq2)[0]

                st.metric(
                    "Global alignment score",
                    f"{alignment.score:.2f}"
                )

                st.code(str(alignment)[:12000])

                st.caption(
                    "This is a global pairwise alignment, not a BLASTp "
                    "database search. The score is not percent identity."
                )

            except ImportError:
                st.error(
                    "Biopython is missing. Add biopython to "
                    "requirements.txt and redeploy."
                )


# ==================================================
# RESULTS & VISUALIZATION
# ==================================================

elif page == "Results & Visualization":
    st.header("📊 Results & Visualization")

    if "blast_results" in st.session_state:
        rows = st.session_state["blast_results"]

        if rows:
            df = pd.DataFrame(rows)

            col1, col2, col3 = st.columns(3)

            col1.metric("Reported hits", len(df))
            col2.metric(
                "Highest identity",
                f"{df['Identity (%)'].max():.2f}%"
            )
            col3.metric(
                "Best E-value",
                f"{df['E-value'].min():.2g}"
            )

            st.dataframe(
                df,
                use_container_width=True,
                hide_index=True
            )

            st.download_button(
                "Export results",
                df.to_csv(index=False),
                file_name="allergenlens_results.csv",
                mime="text/csv"
            )
        else:
            st.info("No BLASTp hits are available.")
    else:
        st.info(
            "Run a BLASTp search to display live result statistics here."
        )


# ==================================================
# EVIDENCE & LIMITATIONS
# ==================================================

elif page == "Evidence & Limitations":
    st.header("📚 Evidence & Limitations")

    st.markdown("""
    **Data sources**
    - UniProtKB for protein records and sequences.
    - NCBI BLASTp for database similarity searches.

    **Reproducibility**
    - Preserve accessions, BLAST request IDs, database names,
      and search dates.
    - Save downloaded results and report search parameters.
    - Search results can change as public databases are updated.

    **Interpretation**
    - Sequence similarity does not establish clinical
      cross-reactivity or individual allergy risk.
    - A database match does not automatically confirm that
      a protein is a clinically established allergen.
    - Pairwise alignment scores and BLASTp scores are different
      and should not be treated as interchangeable.
    """)

st.divider()
st.caption(
    "AllergenLens | Computational Bioinformatics Platform"
)
