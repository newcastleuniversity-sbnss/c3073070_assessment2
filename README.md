# Assignment 2

Technical Report


## Author

c3073070, Newcastle University 2024


## Process 

R script in the repository (R_Script.sh) should be fully executable in R providing the count data also provided in the repository is present in the working directory.

The script outlines the following process:

* Data installation
* Preparation of count data
* Data Quality control (Dispersion estimation, PCA, Heirarchical clustering, heatmap)
* Differential gene expression analysis using DESeq2 package
* Extraction of gene annotations from paired transcriptome database
* Generation of Volcano plots

### Differential Gene Expression Analysis of RNA-seq data

The following analysis is conducted using the DESeq2 package for downsteam differential gene expression analysis of  the publicly available dataset [GSE116583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116583 "@embed") obtained from Gene Expression Omnibus (GEO).

Information about the dataset: https://doi.org/10.1111/ajt.15751


## Data

Download for: 

[Count data](https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/results/counts.zip "@embed") used in script

[RNA-Seq data for GSE116538](https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/results/counts.zip "@embed")




## Results
Images generated for the results section are stored in the Results folder accessible in this repository.


## Software

Salmon (v1.9.0)

R (v4.3.2)
