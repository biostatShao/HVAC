# HVAC
![image](https://github.com/biostatShao/HVAC/blob/main/image.jpg)
Hierarchical Bayesian model with Variational inference and functional Annotation integration for Cross-ancestry prediction
===========================================================================
# 1 Overview

We introduce HVAC  (Hierarchical Bayesian model with Variational inference and functional Annotation integration for Cross-ancestry prediction), a novel cross-population PRS framework employing a three-tiered Bayesian architecture that facilitates coordinated information transfer across diverse populations. The hierarchical structure operates through three functionally integrated components. Firstly, an ancestry-specific inference layer employing annotation-informed Gaussian mixture models to derive population-specific variational posteriors, addressing localized linkage disequilibrium patterns while enabling regularized parameter sharing. Secondly, a cross-ancestry calibration layer utilizing dynamic beta-shrinkage priors that bidirectionally propagate uncertainty estimates across populations, thereby balancing ancestral specificity with shared genetic architecture. Thirdly, polygenic synthesis layer integrating population-specific and cross-ancestry signals through an ensemble super-learner approach with genotype-optimized weighting, employing stacked regression to enhance generalizability while reducing ancestry-related overfitting.

# 2 Installation Instructions
Our software is built on Python version 3.8 or later. To ensure optimal performance and compatibility, we strongly recommend creating a dedicated environment using Conda before installation.
### 2.1 Setting Up the Environment
1. Install Conda: If you do not have Conda installed, download and install it from the official Conda website (https://repo.anaconda.com/archive/index.html).
2. Create a Conda Environment:
```python
conda create -n hvac python=3.11
```
3. Activate the Environment:
```python
conda activate hvac
```
4. Installing Dependencies:
Install `magenpy` using pip
```python
pip install magenpy
```
### 2.2 Building Extensions
1. Ensure the Conda environment is activated.
2. Navigate to the directory containing the setup.py file for the software.
3. Run the following command to build the extensions in place:
```python
python setup.py build_ext --inplace
```

## 3 Usage Instructions

### 3.1 Input Files

#### 3.1.1 Compute Linkage-Disequilibrium (LD) matrices using `magenpy`

The magenpy_ld script is used to compute Linkage-Disequilibrium (LD) matrices, which record the pairwise SNP-by-SNP correlations from a sample of genotype data stored in plink's BED format. The script offers an interface to compute LD matrices by simply specifying the path to the genotype files, the type of LD estimator to use, the subset of variants or samples to keep, and the output directory.
The command processes genotype data for chromosome 1 of the EAS population：
```python
/home/bin/magenpy_ld \
    --bfile ./1KG/sas/hm3_cM_chr1 \
    --min-mac 5 \
    --backend xarray \
    --estimator block \
    --ld-blocks ./EAS_block_coords.bed \
    --storage-dtype int16 \
    --compressor zstd \
    --compression-level 9 \
    --compute-spectral-properties \
    --temp-dir ./temp \
    --output-dir ./SAS/
```

#### 3.1.2 Annotation File

##### 3.1.2.1

First, you need to infer the heritability of each SNP using the S-LDSC model (Baseline-LD annotation can be downloaded [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/)) and use summary statistics that will later be used for training. When running S-LDSC, make sure to use the `--print-coefficients` flag to obtain the regression coefficients.

##### 3.1.2.2

After running S-LDSC:

a. Obtain the `h2g` estimate from the `*.log` file.

b. Retrieve the regression coefficients from the `*.results` file (8th column). Divide the regression coefficients by `h2g` and define it as `T=tau/h2g`. This results in a vector of dimension `C*1`, where `C` is the total number of annotations.

c. From the Baseline-LD annotation downloaded [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/), read the annotation file `baselineLD.*.annot.gz` and retain only the annotation columns (i.e., remove the first 4 columns). This matrix is denoted as `X` with dimensions `M*C`, where `M` is the number of SNPs and `C` is the total number of annotations.

d. Define the expected heritability of each SNP as an `M*1` vector, which is the result of matrix `X` multiplied by `T`.

The final file format is:
- Column 1: SNP ID
- Column 2: Heritability of each SNP

#### 3.1.3 GWAS Summary Statistics File

The file must adhere strictly to the following format:
- `CHR`: Chromosome
- `SNP`: SNP ID
- `BP`: Physical position (base pairs)
- `A1`: Minor allele name
- `A2`: Major allele name
- `P`: GWAS p-value
- `BETA`: Effect size
- `Z`: Z-statistic

#### 3.1.4 Phenotype File

Format:
- Column 1: UDI
- Column 2: Phenotype
- Column 3- (optional): Covariates

### 3.2 Input Parameters

- `traits`: Trait name (e.g., Height). The prefix of the summary data file should be consistent (e.g., Height_sums.txt). Directory name for output files.
- `chr`: Chromosome number (e.g., 1)
- `N`: The GWAS sample size
- `h2`: Estimated SNP Heritability (pre-compute this using your favorite method).
- `sums_p`: Absolute path to the summary data (e.g., `/your/path/sum_file/`)
- `base_p`: Prefix of Plink format LD data files, without the chromosome number (e.g., for chromosome 1, the file is `/your/path/plink/eur_hm3_chr1`, so enter `/your/path/plink/eur_hm3_chr`)
- `target_p`: Prefix of Plink format test set data files, without the chromosome number (e.g., for chromosome 1, the file is `/your/path/plink/ukb22828_eur_hm3_chr1`, so enter `/your/path/plink/ukb22828_eur_hm3_chr`)
- `pheno`: Phenotype file and its absolute path. If covariates are included, they should be in this file (e.g., `/your/path/ukbb.phen`)
- `phe_trait`: Column name of the outcome in the phenotype file (e.g., Height)
- `out`: Output path for result files (e.g., `/your/path/out/`)
- `temp`: Output path for temporary files (e.g., `/your/path/temp/`)
- `cova`: Covariates to consider; if none, enter `NULL` (e.g., `c("BaseAge","Sexgenetic")`)
- `bina`: Whether the outcome is binary data (e.g., `T`)
- `ct_result`: Whether to output 11 tissue-specific PRS prediction levels (e.g., `T`)
- `software_path`: The directory to OmniPRS (e.g., `/your/path/`)
- `plink_path`: The directory to PLINK1.9 and PLINK2.0

## 3.3 Running

Below is an example of running the OmniPRS package:

```r
# R version
R 4.1.2

# Set parameters
base_p <- "/data2/projects/bioinfo/1KG/"
base_f <- "eur_hm3_chr"
target_p <- "/data2/projects/bioinfo/UKB/"
target_f <- "ukb22828_eur_hm3_chr"
pheno <- "/data2/projects/bioinfo/ukbb.phen"
temp <- "/data2/projects/bioinfo/temp_ukb/"
cova <- c("BaseAge", "Sexgenetic", paste0("PC", 1:10))
traits <- "Height"
bina <- FALSE
sums_p <- "/data2/projects/bioinfo/sums/"
out <- "/data2/projects/bioinfo/out/"
phe_trait <- traits
print(traits)

# Set GWAS summary file path
sums_p <- paste0(sums_p, traits, "/")

# Run HVAC for each chromosome
for (chr in 1:22) {
  HVAC(traits, chr, sums_p,
    base_p, base_f, cova, target_p, target_f,
    pheno, phe_trait, out, temp, bina)
}
```

### 3.4 Output Files

#### 3.4.1 Posterior Estimates of Effect Sizes ("OmniPRS_beta")

- **Column 1**: SNP ID
- **Column 2**: Minor allele name
- **Columns 3-13**: Posterior estimates of tissue-specific effect sizes for each tissue

#### 3.4.2 Polygenic Risk Scores for Validation Set Individuals ("OmniPRS_score")

- **Column 1**: Individual ID
- **Columns 2-12**: Tissue-specific polygenic risk scores for each tissue
- **Column 13**: Polygenic risk score using the EW algorithm
- **Column 14**: Polygenic risk score using the LASSO algorithm
- **Column 15**: Polygenic risk score using the BMA algorithm

## Acknowledgments
Special thanks to Shadi[https://github.com/shz9/viprs] for providing the foundational code for this software.

## Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Zhonghe Shao](https://github.com/biostatShao) via zhonghe@hust.edu.cn.

## Update
2024-06-17 OmniPRS version 1.0.0
