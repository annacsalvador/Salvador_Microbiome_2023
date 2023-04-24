# Quantile Normalization via preprocessCore
# Cleaned for upload to Github 4/24/23 by ACS

rm(list=ls())
setwd("/Volumes/Transcend\ 1/qtl_paper2_R/data/raw")
pheno <- read.csv("pheno_ASV_correctedCMM_22.csv", row.names = 1)

#######################################################################################################

# Microbiome data
    # Based on Huda Example:

# install if necessary
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("preprocessCore")

# load package
library(preprocessCore)

# Example
# create a matrix using the same example
#mat <- matrix(c(5,2,0,0,4,1,4,0,0,4,0,0),
#              ncol=3)
#mat
#mat.norm <- normalize.quantiles(mat)

colnames(pheno)


# Phyla (D1)
subset <- pheno[,24:29]
subset.names <- colnames(subset)
res <- c()
res <- as.data.frame(matrix(,nrow(subset),))
for(i in subset.names){
  mat <- as.matrix(subset)
  mat.norm <- as.data.frame(normalize.quantiles(mat))
  colnames(mat.norm) <- subset.names
  rownames(mat.norm) <- rownames(subset)
}

#write.csv(mat.norm, "matrix_cmm_normPhyla.csv")

# Genus (D5)
subset <- pheno[,30:87]
subset.names <- colnames(subset)
res <- c()
res <- as.data.frame(matrix(,nrow(subset),))
for(i in subset.names){
  mat <- as.matrix(subset)
  mat.norm <- as.data.frame(normalize.quantiles(mat))
  colnames(mat.norm) <- subset.names
  rownames(mat.norm) <- rownames(subset)
}

#write.csv(mat.norm, "matrix_cmm_normGenus.csv")

# Species (D6)
subset <- pheno[,88:166]
subset.names <- colnames(subset)
res <- c()
res <- as.data.frame(matrix(,nrow(subset),))
for(i in subset.names){
  mat <- as.matrix(subset)
  mat.norm <- as.data.frame(normalize.quantiles(mat))
  colnames(mat.norm) <- subset.names
  rownames(mat.norm) <- rownames(subset)
}

#write.csv(mat.norm, "matrix_cmm_D6_normSpecies.csv")


# For the ASV Data 
subset <- pheno[,170:ncol(pheno)] # there are MUCH fewer phe columns because of cmm for ASV
colnames(subset)[1:10] # Should start with fbfac7ca3133511d2bc1084a09656368
colnames(subset)[ncol(subset)] # Should end with cbab5f434bd6c0322f7506d3416442eb
subset.names <- colnames(subset)
res <- c()
res <- as.data.frame(matrix(,nrow(subset),))
for(i in subset.names){
  mat <- as.matrix(subset)
  mat.norm <- as.data.frame(normalize.quantiles(mat))
  colnames(mat.norm) <- subset.names
  rownames(mat.norm) <- rownames(subset)
}

#write.csv(mat.norm, "matrix_cmm_normASV.csv")


# merge

phyla <- read.csv("matrix_cmm_normPhyla.csv", row.names = 1)
genus <- read.csv("matrix_cmm_normGenus.csv", row.names = 1)
species <- read.csv("matrix_cmm_D6_normSpecies.csv", row.names = 1)
asv <- read.csv("matrix_cmm_normASV.csv", row.names = 1)

