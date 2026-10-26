
import streamlit as st

st.set_page_config(
    page_title="AllergenLens",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 AllergenLens")
st.caption(
    "Food-Allergen Protein Exploration and Sequence Analysis"
)

st.sidebar.title("Navigation")

page = st.sidebar.radio(
    "Go to",
    [
        "Home",
        "Database Explorer",
        "Sequence Laboratory",
        "Protein Comparison",
        "Results & Visualization",
        "Evidence & Limitations"
    ]
)

if page == "Home":
    st.header("Welcome to AllergenLens")

    st.write(
        "Explore food-allergen proteins, retrieve sequences, "
        "perform similarity searches, compare proteins, "
        "and visualize bioinformatics results."
    )

    col1, col2, col3 = st.columns(3)

    col1.metric("Core Modules", "6")
    col2.metric("Sequence Analysis", "BLASTp")
    col3.metric("Data Source", "UniProt")

    st.subheader("Platform Features")

    st.markdown("""
    - Search and explore protein records
    - Retrieve protein sequences
    - Validate FASTA input
    - Analyze sequence similarity
    - Compare protein sequences
    - Visualize and export results
    """)

elif page == "Database Explorer":
    st.header("Allergen Database Explorer")
    st.info("UniProt database integration is our next step.")

elif page == "Sequence Laboratory":
    st.header("Sequence Analysis Laboratory")
    st.text_area("Paste a protein FASTA sequence")
    st.info("Sequence validation and BLASTp will be added next.")

elif page == "Protein Comparison":
    st.header("Protein Comparison")
    st.info("Pairwise sequence comparison will be implemented here.")

elif page == "Results & Visualization":
    st.header("Results & Visualization")
    st.info("Interactive charts and downloadable results will be added here.")

elif page == "Evidence & Limitations":
    st.header("Evidence & Limitations")

    st.markdown("""
    - Record database accessions and retrieval dates.
    - Document sequence-search parameters.
    - Interpret similarity results carefully.
    - Sequence similarity alone does not establish clinical
      cross-reactivity or individual allergy risk.
    """)

st.divider()
st.caption("AllergenLens | Bioinformatics Research Project")
