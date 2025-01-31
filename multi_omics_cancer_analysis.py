import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import pysam

# Inputs
RNA_SEQ_FILE = "rna_seq_data.csv"  # Gene expression matrix (genes x samples)
DNA_SEQ_FILE = "dna_seq.vcf"       # Variant call format file for DNA-seq
METHYLATION_FILE = "methylation_data.csv"  # Methylation beta values (CpGs x samples)
CLINICAL_DATA = "clinical_data.csv"  # Clinical data (samples x clinical features)

# Outputs
OUTPUT_DIR = "results/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load data
print("Loading data...")
rna_seq = pd.read_csv(RNA_SEQ_FILE, index_col=0)
dna_seq = pysam.VariantFile(DNA_SEQ_FILE)
methylation = pd.read_csv(METHYLATION_FILE, index_col=0)
clinical_data = pd.read_csv(CLINICAL_DATA, index_col=0)

# Ensure samples match across datasets
common_samples = set(rna_seq.columns) & set(methylation.columns) & set(clinical_data.index)
rna_seq = rna_seq[common_samples]
methylation = methylation[common_samples]
clinical_data = clinical_data.loc[common_samples]

# Step 1: Differential Expression Analysis (RNA-seq)
def differential_expression_analysis(rna_seq, clinical_data):
    print("Performing differential expression analysis...")
    groups = clinical_data['group']  # Assume 'group' column has 'normal' and 'cancer'
    cancer_samples = groups == 'cancer'
    normal_samples = groups == 'normal'
    
    results = []
    for gene in rna_seq.index:
        cancer_expr = rna_seq.loc[gene, cancer_samples]
        normal_expr = rna_seq.loc[gene, normal_samples]
        t_stat, p_val = ttest_ind(cancer_expr, normal_expr)
        results.append((gene, t_stat, p_val))
    
    results_df = pd.DataFrame(results, columns=["Gene", "T-statistic", "P-value"])
    results_df["FDR"] = multipletests(results_df["P-value"], method="fdr_bh")[1]
    results_df.to_csv(f"{OUTPUT_DIR}/differential_expression_results.csv", index=False)
    return results_df

de_results = differential_expression_analysis(rna_seq, clinical_data)

# Step 2: Novel NGS Sequencing Analysis (DNA-seq)
def novel_ngs_analysis(dna_seq):
    print("Performing novel NGS sequencing analysis...")
    mutation_counts = {}
    for record in dna_seq:
        gene = record.info.get('GENE', 'unknown')
        if gene not in mutation_counts:
            mutation_counts[gene] = 0
        mutation_counts[gene] += 1
    
    mutation_df = pd.DataFrame(mutation_counts.items(), columns=["Gene", "Mutation_Count"])
    mutation_df.to_csv(f"{OUTPUT_DIR}/mutation_counts.csv", index=False)
    return mutation_df

mutation_df = novel_ngs_analysis(dna_seq)

# Step 3: Methylation Analysis
def methylation_analysis(methylation, clinical_data):
    print("Performing methylation analysis...")
    groups = clinical_data['group']
    cancer_samples = groups == 'cancer'
    normal_samples = groups == 'normal'
    
    results = []
    for cpg in methylation.index:
        cancer_meth = methylation.loc[cpg, cancer_samples]
        normal_meth = methylation.loc[cpg, normal_samples]
        u_stat, p_val = mannwhitneyu(cancer_meth, normal_meth)
        results.append((cpg, u_stat, p_val))
    
    results_df = pd.DataFrame(results, columns=["CpG", "U-statistic", "P-value"])
    results_df["FDR"] = multipletests(results_df["P-value"], method="fdr_bh")[1]
    results_df.to_csv(f"{OUTPUT_DIR}/methylation_results.csv", index=False)
    return results_df

meth_results = methylation_analysis(methylation, clinical_data)

# Step 4: Multi-omics Integration (PCA)
def multi_omics_integration(rna_seq, methylation):
    print("Performing multi-omics integration...")
    combined_data = pd.concat([rna_seq.T, methylation.T], axis=1)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(combined_data)
    
    pca = PCA(n_components=2)
    pca_results = pca.fit_transform(scaled_data)
    pca_df = pd.DataFrame(pca_results, columns=["PC1", "PC2"], index=combined_data.index)
    pca_df.to_csv(f"{OUTPUT_DIR}/pca_results.csv")
    return pca_df

pca_df = multi_omics_integration(rna_seq, methylation)

# Step 5: Visualization
def plot_results(de_results, mutation_df, meth_results, pca_df, clinical_data):
    print("Plotting results...")
    # Differential Expression Volcano Plot
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x="T-statistic", y="-log10(P-value)", data=de_results.assign(**{"-log10(P-value)": -np.log10(de_results["P-value"])}))
    plt.title("Volcano Plot of Differential Expression")
    plt.savefig(f"{OUTPUT_DIR}/volcano_plot.png")
    plt.close()

    # Mutation Count Bar Plot
    plt.figure(figsize=(10, 6))
    sns.barplot(x="Gene", y="Mutation_Count", data=mutation_df.sort_values("Mutation_Count", ascending=False).head(20))
    plt.title("Top 20 Mutated Genes")
    plt.xticks(rotation=90)
    plt.savefig(f"{OUTPUT_DIR}/mutation_barplot.png")
    plt.close()

    # Methylation Boxplot
    plt.figure(figsize=(10, 6))
    sns.boxplot(x="CpG", y="U-statistic", data=meth_results.sort_values("U-statistic", ascending=False).head(20))
    plt.title("Top 20 Differentially Methylated CpGs")
    plt.xticks(rotation=90)
    plt.savefig(f"{OUTPUT_DIR}/methylation_boxplot.png")
    plt.close()

    # PCA Plot
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x="PC1", y="PC2", hue=clinical_data['group'], data=pca_df)
    plt.title("PCA of Multi-omics Data")
    plt.savefig(f"{OUTPUT_DIR}/pca_plot.png")
    plt.close()

plot_results(de_results, mutation_df, meth_results, pca_df, clinical_data)

print("Analysis complete! Results saved to the 'results/' directory.")
