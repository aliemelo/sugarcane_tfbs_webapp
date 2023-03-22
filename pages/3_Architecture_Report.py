import streamlit as st
import subprocess
import tempfile
import os


# Functions ------------------------------------------------ #
if "run_script" not in st.session_state:
    st.session_state.run_script = False

def callback():
    st.session_state.run_script = True

# Main --------------------------------------------------- #
st.markdown("# Architecture Report")
use_tree = st.checkbox("Use dendrogram")
use_repeats = st.checkbox("Consider motif repeats")

# Temp file to store tree file
treedata = st.session_state["treedata"]
tmp_tree = tempfile.NamedTemporaryFile(delete=False)
with open(tmp_tree.name, "w") as fout:
    fout.write(treedata)

# Temp files to store input motif file
my_genes_motifs_df = st.session_state["my_genes_motifs_df"]
tmp_gm = tempfile.NamedTemporaryFile()
my_genes_motifs_df.to_csv(tmp_gm.name, header=False, index=False, sep="\t")

# Run architecture script
if st.button("Run architecture report script", on_click=callback):
    # Filter out genes with too many motifs
    perlscript = 'scripts/calculaNumArqEFiltra.pl'
    p = subprocess.Popen(['perl', perlscript, '1000000', '60000000'], stdin=tmp_gm, stdout=subprocess.PIPE)
    tmp_input = tempfile.NamedTemporaryFile()
    with open(tmp_input.name, "w") as fout:
        fout.write(p.stdout.read().decode('utf-8'))

    # Run architecture script
    perl_script = 'scripts/subsetsWithExpressionTreeConfigurableAndUsingPositionsV3.pl'
    if use_tree and use_repeats:
        p2 = subprocess.Popen(['perl', perl_script, '-r', '-t', tmp_gm.name], stdin=tmp_input, stdout=subprocess.PIPE)
    elif use_tree and not use_repeats:
        p2 = subprocess.Popen(['perl', perl_script, '-t', tmp_gm.name], stdin=tmp_input, stdout=subprocess.PIPE)
    elif use_repeats and not use_tree:
        p2 = subprocess.Popen(['perl', perl_script, '-r'], stdin=tmp_input, stdout=subprocess.PIPE)
    else:
        p2 = subprocess.Popen(['perl', perl_script], stdin=tmp_input, stdout=subprocess.PIPE)

    st.sidebar.download_button(label="Download Output", data=p2.stdout.read(), mime="text/plain", file_name="out.txt", on_click=callback)