# OneK1K dataset Analysis
This repository contains the results of my internship at the Division of Computational Genomics and Systems Genetics (DKFZ, Prof. Dr. Stelge, supervisor: Dr. Ueltzhöffer).

Contents:
- `Explorative-analysis.ipynb`: plots comparing RNA expression of single cells with different eQTLs. The eQTLs were selected either by q-value from the eQTL analysis or by the results of the t-test.
- `Wasserstein-distance_KS-test.ipynb`: plots comparing RNA expression of single cells with different eQTLs. The eQTLs were selected either by Wasserstein distance or by the results of the Kolmogorov-Smirnov test.
- `Metacells_CD4NC.ipynb`: plots comparing results of RNA expression of metacells. Moreover, results with different sizes of the metacells were compared.
- `Cumulative_effect_pro_gene.ipynb`: plots comparing cumulative effect sizes of eQTLs acting on the same gene.
- `Gauss-mixture-model.ipynb`: results of modeling the distribution of single-cell RNA expression values with a Gaussian mixture model.
- `ZINB-model.ipynb`: results of modeling the distribution of single-cell RNA expression values with a Zero-Inflated Negative Binomial model.
- `utils_data.py`: script for faster loading of the data for the `Gauss-mixture-model.ipynb` and `ZINB-model.ipynb` notebooks.
