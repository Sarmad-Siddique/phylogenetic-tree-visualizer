import streamlit as st

st.title("Interactive Phylogenetic Tree Builder")

uploaded_file = st.file_uploader("Upload a FASTA file")

if uploaded_file:
    st.success("File uploaded. Proceed to alignment and tree generation.")