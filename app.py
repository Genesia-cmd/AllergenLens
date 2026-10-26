
import streamlit as st
import pandas as pd
import requests
import re
from collections import Counter

st.set_page_config(
    page_title="AllergenLens",
    page_icon="🧬",
    layout="wide"
)

API_URL = "https://rest.uniprot.org/uniprotkb/search"

st.title("🧬 AllergenLens")
st.caption(
    "Interactive Food-Allergen Protein Exploration "
    "and Sequence Analysis Platform"
)

# ---------------- SIDEBAR ----------------
page = st.sidebar.radio(
    "Navigation",
    [
        "Home",
        "Database Explorer",
        "Sequence Laboratory",
        "Protein Comparison",
        "Results & Visualization",
        "Evidence & Limitations"
    ]
)


# ---------------- HELPER FUNCTIONS ----------------
@st.cache_data(ttl=3600)
def search_uniprot(query, limit):
    params = {
        "query": query,
        "format": "json",
        "size": limit,
        "fields": (
            "accession,protein_name,organism_name,"
            "gene_names,length,sequence"
        )
    }

    response = requests.get(
        API_URL,
        params=params,
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


# ---------------- HOME ----------------
if page == "Home":
    st.header("Welcome to AllergenLens")

    st.write(
        "Explore food-allergen proteins, retrieve biological "
        "sequences, investigate sequence similarity, compare "
        "proteins, and visualize computational analysis results."
    )

    col1, col2, col3 = st.columns(3)

    col1.metric("Core Modules", "6")
    col2.metric("Protein Database", "UniProt")
    col3.metric("Similarity Search", "BLASTp")

    st.subheader("Platform Features")

    st.markdown("""
    - Search protein records from public databases
    - Retrieve protein sequences in FASTA format
    - Validate protein sequences
    - Perform protein similarity searches
    - Compare protein sequences
    - Visualize and export analysis results
    """)

    st.info(
        "Use the sidebar to explore the available modules."
    )


# ---------------- DATABASE EXPLORER ----------------
elif page == "Database Explorer":
    st.header("🔎 Allergen Database Explorer")

    st.write(
        "Search UniProtKB for protein records by name, "
        "organism, or accession."
    )

    st.caption(
        "Search matches may include proteins that are not "
        "confirmed food allergens."
    )

    search_term = st.text_input(
        "Protein, allergen, organism, or accession",
        placeholder="e.g. Ara h 2, peanut, P02754"
    )

    limit = st.selectbox(
        "Maximum records",
        [5, 10, 20, 50],
        index=1
    )

    if st.button("Search UniProt", type="primary"):
        if not search_term.strip():
            st.warning("Please enter a search term.")
        else:
            with st.spinner("Searching UniProtKB..."):
                try:
                    records = search_uniprot(
                        search_term.strip(),
                        limit
                    )

                    rows = []

                    for record in records:
                        accession = record.get(
                            "primaryAccession", ""
                        )

                        organism = record.get(
                            "organism", {}
                        ).get("scientificName", "Unknown")

                        sequence_data = record.get("sequence", {})

                        rows.append({
                            "Accession": accession,
                            "Protein": get_protein_name(record),
                            "Organism": organism,
                            "Length (aa)": sequence_data.get("length"),
                            "Reviewed": record.get(
                                "entryType", ""
                            ) == "UniProtKB reviewed (Swiss-Prot)"
                        })

                    if rows:
                        st.session_state["uniprot_records"] = records
                        st.session_state["uniprot_table"] = pd.DataFrame(
                            rows
                        )
                        st.success(
                            f"Retrieved {len(rows)} records."
                        )
                    else:
                        st.session_state.pop(
                            "uniprot_records", None
                        )
                        st.session_state.pop(
                            "uniprot_table", None
                        )
                        st.info("No matching records found.")

                except requests.RequestException as error:
                    st.error(f"UniProt request failed: {error}")
                except (ValueError, KeyError, TypeError) as error:
                    st.error(f"Could not process results: {error}")

    if "uniprot_table" in st.session_state:
        df = st.session_state["uniprot_table"]
        records = st.session_state["uniprot_records"]

        st.subheader("Search Results")
        st.dataframe(
            df,
            use_container_width=True,
            hide_index=True
        )

        st.download_button(
            "Download results as CSV",
            df.to_csv(index=False),
            file_name="allergenlens_proteins.csv",
            mime="text/csv"
        )

        st.subheader("Protein Details")

        options = {
            f"{r.get('primaryAccession', '')} — "
            f"{get_protein_name(r)}": i
            for i, r in enumerate(records)
        }

        selected = st.selectbox(
            "Select a protein",
            list(options.keys())
        )

        record = records[options[selected]]

        accession = record.get("primaryAccession", "")
        protein_name = get_protein_name(record)

        organism = record.get(
            "organism", {}
        ).get("scientificName", "Unknown")

        sequence = record.get("sequence", {}).get("value", "")

        st.markdown(f"**Protein:** {protein_name}")
        st.markdown(f"**Accession:** `{accession}`")
        st.markdown(f"**Organism:** {organism}")
        st.markdown(f"**Sequence length:** {len(sequence)} amino acids")

        if sequence:
            fasta = f">{accession} {protein_name}\n"
            fasta += "\n".join(
                sequence[i:i + 60]
                for i in range(0, len(sequence), 60)
            )

            with st.expander("View FASTA sequence"):
                st.code(fasta, language=None)

            st.download_button(
                "Download FASTA",
                fasta,
                file_name=f"{accession}.fasta",
                mime="text/plain",
                key="download_selected_fasta"
            )

        st.markdown(
            f"[View record on UniProt]"
            f"(https://www.uniprot.org/uniprotkb/{accession}/entry)"
        )


# ---------------- SEQUENCE LABORATORY ----------------
elif page == "Sequence Laboratory":
    st.header("🧪 Sequence Analysis Laboratory")

    st.write(
        "Upload a protein FASTA file or paste an amino-acid "
        "sequence to validate it before similarity analysis."
    )

    uploaded_file = st.file_uploader(
        "Upload protein FASTA",
        type=["fasta", "fa", "faa", "txt"],
        key="protein_upload"
    )

    pasted_sequence = st.text_area(
        "Or paste a protein sequence / FASTA",
        placeholder=">Example_protein\nMKWVTFISLL...",
        height=160,
        key="protein_input"
    )

    if uploaded_file is not None:
        raw_text = uploaded_file.getvalue().decode(
            "utf-8", errors="replace"
        )
        source = "Uploaded file"
    else:
        raw_text = pasted_sequence
        source = "Pasted sequence"

    if st.button("Validate Protein Sequence", type="primary"):
        lines = [
            line.strip()
            for line in raw_text.splitlines()
            if line.strip()
        ]

        headers = [
            line for line in lines
            if line.startswith(">")
        ]

        if not lines:
            st.error("Please upload or paste a sequence.")

        elif len(headers) > 1:
            st.error(
                "Please provide one protein sequence at a time."
            )

        else:
            sequence_lines = [
                line for line in lines
                if not line.startswith(">")
            ]

            sequence = re.sub(
                r"\s+", "", "".join(sequence_lines)
            ).upper()

            allowed = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")

            if not sequence:
                st.error("No amino-acid sequence was found.")

            else:
                sequence = sequence.rstrip("*")

                invalid = sorted(set(sequence) - allowed)

                if invalid:
                    st.error(
                        "Invalid protein symbols: "
                        + ", ".join(invalid)
                    )

                elif len(sequence) < 10:
                    st.warning(
                        "The sequence is very short. "
                        "Check whether it is suitable for your analysis."
                    )

                else:
                    st.session_state["validated_sequence"] = sequence

                    st.success(
                        "Protein sequence validation passed."
                    )

                    col1, col2, col3 = st.columns(3)

                    col1.metric(
                        "Sequence Length",
                        f"{len(sequence):,} aa"
                    )

                    col2.metric(
                        "Unique Symbols",
                        len(set(sequence))
                    )

                    col3.metric("Input Source", source)

                    counts = Counter(sequence)

                    composition = pd.DataFrame({
                        "Amino Acid": list(counts.keys()),
                        "Count": list(counts.values())
                    }).sort_values("Amino Acid")

                    st.subheader("Amino-Acid Composition")

                    st.bar_chart(
                        composition.set_index("Amino Acid")["Count"]
                    )

                    ambiguous = set("BXZJUO")
                    ambiguous_count = sum(
                        counts[aa] for aa in ambiguous
                    )

                    if ambiguous_count:
                        st.warning(
                            f"{ambiguous_count} residue(s) contain "
                            "ambiguous or uncommon amino-acid symbols."
                        )

                    header = (
                        headers[0]
                        if headers
                        else ">AllergenLens_query"
                    )

                    fasta_output = header + "\n" + "\n".join(
                        sequence[i:i + 60]
                        for i in range(0, len(sequence), 60)
                    )

                    with st.expander("View Validated FASTA"):
                        st.code(fasta_output, language=None)

                    st.download_button(
                        "Download Validated FASTA",
                        fasta_output,
                        file_name="allergenlens_query.fasta",
                        mime="text/plain"
                    )

                    st.info(
                        "Sequence validation is complete. "
                        "It has not yet been submitted to BLASTp."
                    )

    if "validated_sequence" in st.session_state:
        st.caption(
            "A validated sequence is available in this session."
        )


# ---------------- PROTEIN COMPARISON ----------------
elif page == "Protein Comparison":
    st.header("🧬 Protein Comparison")

    st.write(
        "This module will compare two protein sequences "
        "using alignment-based statistics."
    )

    st.info(
        "Pairwise alignment, percentage identity, query coverage, "
        "and E-value reporting are the next development steps."
    )


# ---------------- RESULTS ----------------
elif page == "Results & Visualization":
    st.header("📊 Results & Visualization")

    st.write(
        "This module will present protein search results, "
        "interactive charts, and downloadable reports."
    )

    st.info(
        "BLASTp result visualization will be connected "
        "after the live search module is implemented."
    )


# ---------------- EVIDENCE ----------------
elif page == "Evidence & Limitations":
    st.header("📚 Evidence & Limitations")

    st.markdown("""
    **Database source**
    - UniProtKB provides protein records and sequences.
    - Accession IDs allow users to locate the source records.

    **Reproducibility**
    - Record retrieval dates and search parameters.
    - Distinguish stored results from newly calculated results.

    **Scientific limitations**
    - A database search match does not automatically confirm
      that a protein is an allergen.
    - Sequence similarity alone cannot establish clinical
      cross-reactivity or individual allergy risk.
    """)

st.divider()
st.caption(
    "AllergenLens | Computational Bioinformatics Platform"
)
