import streamlit as st
import matplotlib.pyplot as plt
import networkx as nx
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import os

# Step 1: Function to save user sequences to a FASTA file
def save_fasta(user_sequences):
    fasta_file = "sequences.fasta"
    with open(fasta_file, "w") as f:
        f.write(user_sequences)
    return fasta_file

# Step 2: Build a phylogenetic tree
def build_phylogenetic_tree(file_path):
    try:
        alignment = AlignIO.read(file_path, "fasta")
        calculator = DistanceCalculator("identity")
        constructor = DistanceTreeConstructor(calculator, "upgma")
        phylo_tree = constructor.build_tree(alignment)

        # Save the tree as an image
        tree_image_file = "tree.png"
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        Phylo.draw(phylo_tree, axes=ax)
        plt.savefig(tree_image_file)  # Save image
        plt.close()  # Close figure

        return phylo_tree, tree_image_file
    except Exception as e:
        return None, str(e)

# Step 3: Visualize SNP interaction network
def visualize_interaction_network(edges):
    G = nx.Graph()
    G.add_edges_from(edges)
    pos = nx.spring_layout(G)
    fig, ax = plt.subplots()
    nx.draw(G, pos, with_labels=True, node_color="lightblue", edge_color="gray", node_size=2000, ax=ax)
    return fig

# Streamlit UI
st.title("Phylogenetic Tree & Sequence Analysis")

# User input for FASTA sequences
st.subheader("Enter your sequences in FASTA format:")
user_input = st.text_area("Example format:\n"
                          ">Seq1\nATGCGTACGTTAGTAACTG\n"
                          ">Seq2\nATGCGTACGTTAGTACCTG\n"
                          ">Seq3\nATGCGTACGTTGGTAACTG\n"
                          ">Seq4\nATGCGGACGTTAGTAACTG", height=200)

if st.button("Generate Phylogenetic Tree & Analyze Sequences"):
    if user_input.strip() == "":
        st.error("Please enter valid sequences in FASTA format.")
    else:
        fasta_file = save_fasta(user_input)
        st.success("Phylogenetic Tree Generated Successfully!")

        # Build tree
        tree, tree_image = build_phylogenetic_tree(fasta_file)
        if tree and os.path.exists(tree_image):
            st.image(tree_image, caption="Phylogenetic Tree")
        else:
            st.error(f"Error generating tree: {tree_image}")

        # Visualize SNP interaction network
        st.subheader("SNP Interaction Network")
        example_edges = [("Seq1", "Seq2"), ("Seq2", "Seq3"), ("Seq3", "Seq4")]
        fig_network = visualize_interaction_network(example_edges)
        st.pyplot(fig_network)