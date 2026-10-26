
import streamlit as st
import pandas as pd
import requests

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


if page == "Home":
    st.header("Welcome to AllergenLens")

    st.write(
        "Explore food-allergen proteins, retrieve biological "
        "sequences, investigate sequence similarity, and "
        "visualize computational analysis results."
    )

    c1, c2, c3 = st.columns(3)
    c1.metric("Core Modules", "6")
    c2.metric("Protein Database", "UniProt")
    c3.metric("Analysis Method", "BLASTp")

    st.subheader("Platform capabilities")

    st.markdown("""
    - Search protein records from a public database
    - Retrieve protein sequences in FASTA format
    - Validate user-submitted protein sequences
    - Perform protein similarity searches
    - Compare protein sequences and visualize results
    - Export data for further analysis
    """)

elif page == "Database Explorer":
    st.header("🔎 Allergen Database Explorer")

    st.write(
        "Search UniProtKB for protein records. "
        "Search results may include proteins that are not "
        "confirmed allergens."
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
                        sequence = record.get("sequence", {})

                        rows.append({
                            "Accession": accession,
                            "Protein": get_protein_name(record),
                            "Organism": organism,
                            "Length (aa)": sequence.get("length"),
                            "Reviewed": record.get(
                                "entryType", ""
                            ) == "UniProtKB reviewed (Swiss-Prot)"
                        })

                    if rows:
                        st.session_state["uniprot_records"] = records
                        st.session_state["uniprot_table"] = pd.DataFrame(rows)
                        st.success(
                            f"Found {len(rows)} protein records."
                        )
                    else:
                        st.session_state.pop("uniprot_records", None)
                        st.session_state.pop("uniprot_table", None)
                        st.info("No matching records found.")

                except requests.RequestException as error:
                    st.error(f"UniProt request failed: {error}")
                except (ValueError, KeyError, TypeError) as error:
                    st.error(f"Could not process results: {error}")

    if "uniprot_table" in st.session_state:
        df = st.session_state["uniprot_table"]
        records = st.session_state["uniprot_records"]

        st.subheader("Search results")
        st.dataframe(df, use_container_width=True, hide_index=True)

        st.download_button(
            "Download results as CSV",
            df.to_csv(index=False),
            file_name="allergenlens_proteins.csv",
            mime="text/csv"
        )

        st.subheader("Explore an individual protein")

        options = {
            f"{r.get('primaryAccession', '')} — "
            f"{get_protein_name(r)}": i
            for i, r in enumerate(records)
        }

        selected = st.selectbox(
            "Select a protein record",
            list(options.keys())
        )

        record = records[options[selected]]
        sequence = record.get("sequence", {}).get("value", "")
        accession = record.get("primaryAccession", "")
        organism = record.get("organism", {}).get(
            "scientificName", "Unknown"
        )

        st.markdown(f"**Protein:** {get_protein_name(record)}")
        st.markdown(f"**Accession:** `{accession}`")
        st.markdown(f"**Organism:** {organism}")
        st.markdown(f"**Sequence length:** {len(sequence)} amino acids")

        if sequence:
            fasta = f">{accession} {get_protein_name(record)}\n"
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
                mime="text/plain"
            )

        st.markdown(
            f"[View record on UniProt](https://www.uniprot.org/uniprotkb/{accession}/entry)"
        )

elif page == "Sequence Laboratory":
    st.header("🧪 Sequence Laboratory")
    st.info(
        "Next module: FASTA upload, sequence validation, "
        "and live BLASTp submission."
    )

elif page == "Protein Comparison":
    st.header("🧬 Protein Comparison")
    st.info(
        "Next module: pairwise sequence identity, coverage, "
        "E-value, and alignment display."
    )

elif page == "Results & Visualization":
    st.header("📊 Results & Visualization")
    st.info(
        "Next module: interactive charts, downloadable "
        "analysis tables, and validation results."
    )

elif page == "Evidence & Limitations":
    st.header("📚 Evidence & Limitations")

    st.markdown("""
    - Source: UniProtKB protein database.
    - Record accessions should be retained for reproducibility.
    - Record retrieval dates and search settings.
    - Search matches do not automatically establish allergen status.
    - Sequence similarity alone cannot confirm clinical
      cross-reactivity or individual allergy risk.
    """)

st.divider()
st.caption("AllergenLens | Computational Bioinformatics Platform")
