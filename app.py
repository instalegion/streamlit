
import streamlit as st
import scanpy as sc
import os
import matplotlib.pyplot as plt
import io
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from pandas.api.types import CategoricalDtype

st.set_page_config(page_title='scRNA-seq Analysis App', layout='wide')

st.sidebar.title("Επιλογές Ανάλυσης")
tab = st.sidebar.radio("Επιλέξτε λειτουργία:", [
    "Εισαγωγή Δεδομένων",
    "Προεπεξεργασία",
    "DEG & Plots",
    "Πληροφορίες Ομάδας"
])

if tab == "Εισαγωγή Δεδομένων":
    st.title("Εισαγωγή Αρχείων .h5ad")
    uploaded_file = st.file_uploader("Ανεβάστε αρχείο .h5ad", type="h5ad")
    if uploaded_file:
        os.makedirs("data", exist_ok=True)
        input_path = os.path.join("data", "input_data.h5ad")
        with open(input_path, "wb") as f:
            f.write(uploaded_file.read())
        adata = sc.read_h5ad(input_path)
        st.session_state["adata"] = adata
        st.success("Ανέβηκε επιτυχώς")
        st.write(f"{adata.n_obs} κύτταρα - {adata.n_vars} γονίδια")
        st.dataframe(adata.obs.head())

elif tab == "Προεπεξεργασία":
    st.title("Προεπεξεργασία")

    def slider_with_input(label, min_val, max_val, default, step=1):
        col1, col2 = st.columns([2, 1])
        with col1:
            value = st.slider(label, min_val, max_val, default, step)
        with col2:
            value = st.number_input(" ", min_value=min_val, max_value=max_val, value=value, step=step, key=label)
        return value

    if "adata" not in st.session_state:
        st.warning("Ανεβάστε πρώτα αρχείο.")
        st.stop()

    adata = st.session_state["adata"].copy()

    min_genes = slider_with_input("min_genes", 0, 1000, 200)
    min_cells = slider_with_input("min_cells", 1, 20, 3)

    if st.button("Εκκίνηση Προεπεξεργασίας"):
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        adata.raw = adata
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        adata = adata[:, adata.var.highly_variable].copy()
        sc.pp.scale(adata)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata)
        sc.tl.umap(adata)
        st.session_state["adata"] = adata
        st.success("Η προεπεξεργασία ολοκληρώθηκε.")


elif tab == "DEG & Plots":
    st.title("Διαφορική Έκφραση & Γραφήματα")
    if "adata" not in st.session_state:
        st.warning("Χρειάζεται ανεβασμένο και προεπεξεργασμένο αρχείο.")
        st.stop()

    adata = st.session_state["adata"].copy()

    groupby = st.selectbox("Ομαδοποίηση:", [col for col in adata.obs.columns if isinstance(adata.obs[col].dtype, (CategoricalDtype, object.__class__))])
    reference = st.selectbox("Ομάδα Αναφοράς:", list(adata.obs[groupby].dropna().unique()))
    method = st.selectbox("Μέθοδος:", ["wilcoxon", "t-test", "logreg"])

    if st.button("Εκκίνηση DEG"):
        try:
            adata_filtered = adata[adata.obs[groupby].notna()].copy()
            sc.tl.rank_genes_groups(adata_filtered, groupby=groupby, reference=reference, method=method)
            st.success("Η ανάλυση ολοκληρώθηκε.")
            plot_type = st.radio("Επιλογή Γραφήματος:", ["Ranked Genes", "Heatmap", "Dotplot", "Violin", "Volcano"])

            fig = None
            if plot_type == "Ranked Genes":
                fig = sc.pl.rank_genes_groups(adata_filtered, n_genes=10, show=False, return_fig=True)
            elif plot_type == "Heatmap":
                fig = sc.pl.rank_genes_groups_heatmap(adata_filtered, n_genes=10, show=False, return_fig=True)
            elif plot_type == "Dotplot":
                fig = sc.pl.rank_genes_groups_dotplot(adata_filtered, n_genes=10, show=False, return_fig=True)
            elif plot_type == "Violin":
                fig = sc.pl.rank_genes_groups_violin(adata_filtered, n_genes=10, show=False, return_fig=True)
            elif plot_type == "Volcano":
                result = adata_filtered.uns["rank_genes_groups"]
                group = result["names"].dtype.names[0]
                lfc = result["logfoldchanges"][group]
                pvals = result["pvals_adj"][group]
                fig = plt.figure(figsize=(6, 4))
                plt.scatter(lfc, -np.log10(pvals), alpha=0.6)
                plt.xlabel("logFC")
                plt.ylabel("-log10(p-adj)")
                plt.title("Volcano Plot")

            
            if fig:
                if isinstance(fig, list):
                    fig = fig[0]
                if hasattr(fig, 'figure'):
                    fig = fig.figure
                st.pyplot(fig)
                buf = io.BytesIO()
                fig.savefig(buf, format="png", dpi=300)
                buf.seek(0)
                st.download_button("Κατέβασε PNG", data=buf, file_name=f"{plot_type}.png", mime="image/png")

        except Exception as e:
            st.error(f"Σφάλμα: {str(e)}")

elif tab == "Πληροφορίες Ομάδας":
    st.title("Συντελεστές")
    st.markdown("""
    **Ομάδα Έργου scRNA-seq App**  
    - **Σουλιώτης Κωνσταντίνος**: Streamlit UI, προεπεξεργασία  
    - **Αντώνης Σοφικίτης**: DEG ανάλυση, τεκμηρίωση, δομή εφαρμογής  
    """)
