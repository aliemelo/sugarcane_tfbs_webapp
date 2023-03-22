import streamlit as st
import pandas as pd

st.set_page_config(layout="wide")

# Add a title and intro text
st.title("TFBS Analysis in sugarcane SP80-3280")
st.markdown("This is a web app to allow exploration of TFBSs in the sugarcane variety SP80-3280.")

st.markdown("### How to use this web app?")
st.markdown("1. Provide a list of target genes in the *Gene Selection* page")
st.markdown("2. Select the experiments and tissues you want to use in the *Experiment Selection* page")
st.markdown("3. Compute the possible promoter architectures for the target genes in the *Architecture Report* page")
st.markdown("4. Visualise the heatmap for the selected experiments and the distribution of motifs in the promoter region of the target genes in the *Visualisation* page")

