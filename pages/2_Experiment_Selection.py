import streamlit as st
import pandas as pd
from scipy.cluster.hierarchy import linkage, to_tree, dendrogram
from scipy.spatial import distance

# Functions ------------------------------------------------ #
@st.cache_data
def convert_df(df):
    return df.to_csv().encode('utf-8')

def cluster_by_exp(df):
    row_linkage = linkage(distance.pdist(df), method='complete', metric='euclidean')
    dendro_exp = dendrogram(row_linkage, orientation='left', labels=df.index, no_plot=True)
    return row_linkage, dendro_exp

def write2newick(exp_df_genes, row_linkage):
    # write dendrogram in Newick format
    indexes = exp_df_genes.index.tolist()
    tree_exp = to_tree(row_linkage, False)
    treedata = getNewick(tree_exp, "", tree_exp.dist, indexes)
    return treedata

def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick


# Main --------------------------------------------------- #
# Load expression data for the selected genes
conn = st.session_state['conn']
genes = st.session_state['genes']
exp_df = pd.read_sql(f"SELECT * FROM expression_renamed WHERE gene in {genes}", conn)

# Select experimental data
st.markdown("# Experiment selection")

experiments = ['gene']

col1, col2, col3 = st.columns(3)

with col1:
    st.markdown("**Experiment**")
    exp_drought = st.checkbox("Drought")
    exp_ancestral = st.checkbox("Ancestral")
    exp_biomass = st.checkbox("Biomass")

if exp_ancestral or exp_biomass or exp_drought:
    with col2:
        st.markdown("**Tissue**")
        tissue_leaf = st.checkbox("Leaf +1")
        if exp_biomass or exp_ancestral:
            tissue_inter1 = st.checkbox("Internode 1")
            tissue_inter5 = st.checkbox("Internode 5")
        if exp_ancestral:
            tissue_inter9 = st.checkbox("Internode 9")
        if exp_drought:
            tissue_root = st.checkbox("Root")

if exp_biomass or exp_drought:
    with col3:
        if exp_biomass:
            st.markdown("**Collection time**")
            if tissue_leaf or tissue_inter1:
                time_4m = st.checkbox("4 months")
            time_8m = st.checkbox("8 months")
            time_12m = st.checkbox("12 months")
        if exp_drought:
            st.markdown("**Condition**")
            control = st.checkbox("Control")
            treated = st.checkbox("Treated")

if exp_ancestral:
    if tissue_leaf:
        experiments.append('L1_Ancestral')
    if tissue_inter1:
        experiments.append('I1_Ancestral')
    if tissue_inter5:
        experiments.append('I5_Ancestral')
    if tissue_inter9:
        experiments.append('I9_Ancestral')

if exp_biomass:
    if tissue_leaf:
        if time_4m:
            experiments.append('L1_4M_Biomass')
        if time_8m:
            experiments.append('L1_8M_Biomass')
        if time_12m:
            experiments.append('L1_12M_Biomass')
    if tissue_inter1:
        if time_4m:
            experiments.append('I1_4M_Biomass')
        if time_8m:
            experiments.append('I1_8M_Biomass')
        if time_12m:
            experiments.append('I1_12M_Biomass')
    if tissue_inter5:
        if time_8m:
            experiments.append('I5_8M_Biomass')
        if time_12m:
            experiments.append('I5_12M_Biomass')

if exp_drought:
    if tissue_leaf:
        if control:
            experiments.append("L1_Drought_control")
        if treated:
            experiments.append("L1_Drought_treated")
    if tissue_root:
        if control:
            experiments.append("RT_Drought_control")
        if treated:
            experiments.append("RT_Drought_treated")

if not experiments:
    st.error("Please select at least one experiment")
else:
    data = exp_df[experiments]
    data = data.set_index('gene')
    st.session_state["data"] = data
    csv = convert_df(data)
    st.table(data)

    # Option to download the selected expression data
    st.sidebar.markdown("**Download Output**")
    st.sidebar.download_button(label="Download data as CSV", data=csv, mime="text/csv", file_name="genes_exp.csv")
    
    # Calculate dendrogram
    row_linkage, dendro_exp = cluster_by_exp(data)
    treedata = write2newick(data, row_linkage)
    st.session_state["treedata"] = treedata
    st.sidebar.download_button(label="Download dendrogram in Newick format", 
                               data=treedata, mime="text/plain", file_name="tree.newick")