# Multi-Omics-Cancer-Analysis
Multi-omics analysis in cancer research, including deep statistical analysis, novel NGS sequencing analysis, and visualization. The script assumes you have RNA-seq, DNA-seq, and methylation data.

# Multi-Omics Cancer Analysis

This repository contains a Python script for performing multi-omics analysis in cancer research. The script integrates RNA-seq, DNA-seq, and methylation data to identify differentially expressed genes, mutated genes, and differentially methylated CpGs. It also performs multi-omics integration using PCA and visualizes the results.

## Inputs
1. `rna_seq_data.csv`: Gene expression matrix (genes x samples).
2. `dna_seq.vcf`: Variant call format file for DNA-seq.
3. `methylation_data.csv`: Methylation beta values (CpGs x samples).
4. `clinical_data.csv`: Clinical data (samples x clinical features).

## Outputs
1. `differential_expression_results.csv`: Differential expression analysis results.
2. `mutation_counts.csv`: Mutation counts per gene.
3. `methylation_results.csv`: Differential methylation analysis results.
4. `pca_results.csv`: PCA results for multi-omics integration.
5. Visualizations: Volcano plot, mutation bar plot, methylation boxplot, and PCA plot.

## Requirements
- Python 3.8+
- Libraries listed in `requirements.txt`

## Usage
1. Clone the repository.
2. Install dependencies: `pip install -r requirements.txt`.
3. Place input files in the root directory.
4. Run the script: `python multi_omics_cancer_analysis.py`.
5. Check the `results/` directory for outputs.

## License
MIT
