import streamlit as st
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from io import StringIO
import os
import tempfile
from ete3 import Tree, TreeStyle, TextFace
import re

st.set_page_config(page_title="Phylogenetic Tree Builder", layout="centered")
st.title("🔬 Interactive Phylogenetic Tree Builder")

uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta", "fa"])

method = st.selectbox("Choose Tree Construction Method", ["Neighbor-Joining", "UPGMA"])

if uploaded_file is not None:
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "input.fasta")
        aln_path = os.path.join(tmpdir, "aligned.aln")
        dnd_path = os.path.join(tmpdir, "tree.nwk")

        # Save uploaded file
        with open(fasta_path, "wb") as f:
            f.write(uploaded_file.read())

        # Run Clustal Omega with tree output
        clustalomega_cline = ClustalOmegaCommandline(infile=fasta_path, outfile=aln_path, verbose=True, auto=True, force=True, guidetree_out=dnd_path)
        stdout, stderr = clustalomega_cline()

        # Read the alignment
        alignment = AlignIO.read(aln_path, "fasta")
        st.subheader("🧬 Aligned Sequences")
        st.text(str(alignment))

        # Calculate distance matrix
        calculator = DistanceCalculator("identity")  # Replace with Kimura if needed
        distance_matrix = calculator.get_distance(alignment)

        # Build the tree
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(distance_matrix) if method == "Neighbor-Joining" else constructor.upgma(distance_matrix)

        # Convert tree to Newick and clean it for ete3
        newick_io = StringIO()
        Phylo.write(tree, newick_io, "newick")
        newick_str = newick_io.getvalue()

        # Remove inner node names like 'Inner1'
        newick_str = re.sub(r'Inner\d+', '', newick_str)

        # Display interactive tree
        ete_tree = Tree(newick_str)

        def layout(node):
            if node.is_leaf():
                name_face = TextFace(node.name, fsize=10)
                node.add_face(name_face, column=0)
                
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.scale = 120
        ts.layout_fn = layout

        st.subheader("🌳 Interactive Tree (Zoom + Pan)")
        html_content, _ = ete_tree.render("%%return", tree_style=ts)
        st.components.v1.html(html_content, height=600, scrolling=True)
