Multivariate genome-wide association analysis

Multi-variate GWAS (mvGWAS) is an extension of traditional GWAS that analyzes multiple traits simultaneously rather than one trait at a time. This approach increases statistical power and helps detect pleiotropic effects, where a single genetic variant influences multiple traits.

As input, I have plink files (.bim, .bed, .map, .fam, .ped) and phenotype file (mv_gwas_phenotypes.txt)

1) Running Multivariate GWAS in PLINK 2.0
PLINK 2.0 can perform multivariate regression using the --glm multivariate option.
plink2 --bfile alz8 --pheno mv_gwas_phenotypes.txt --glm multivariate --covar mv_gwas_phenotypes.txt --out mv_gwas_results


2) Run Multivariate GWAS Using GEMMA
Since GEMMA supports multivariate linear mixed models (mvLMM), we need to correctly prepare input files before running the analysis.

Before running multivariate GWAS, you must estimate the kinship matrix:

gemma -bfile alz8 -gk 1 -o kinship

This generates:
output/kinship.cXX.txt → Genetic relationship matrix (GRM).
Once the kinship matrix is ready, run GEMMA multivariate GWAS:
gemma -bfile alz8 -p mv_gwas_phenotypes.txt -k output/kinship.cXX.txt -lmm 4 -o mv_gwas_results


