# Core Measureable Microbiota
# Cleaned for upload to Github 4/24/23 by ACS

rm(list = ls())
setwd("/Volumes/Transcend\ 1/qtl_paper2_R/data/raw")

  # Defined as only microbial traits present in 20% of study population individuals used for QTL mapping. 
  # Kemis 2019: Genetic determinants of gut microbiota composition and bile acid profiles in mice
  # They did this step after the 16S analysis so it should happen before I do the quantile normalization with the data

phe <- read.csv(file = "pheno_D1_to_D6.csv", row.names = 1) # contains microbial and metabolic traits from our studies

# only need microbiome data and the main identifiers
colnames(phe)

subset <- phe[,c(1:3,25:ncol(phe))]

# The loop needs to find non-zero data and select columns were at least 20% of the data is non-zero

  # Example of finding zeros:
#example = matrix(sample(c(0,0,0,100),size=70,replace=T),ncol=7)
#example
#colSums(example != 0)

  # Test finding zeros
colsumm <- colSums(subset != 0, na.rm = T) # You have to remove NAs for colSums to work
percentage <- colsumm/nrow(subset)
colsumm20 <- percentage >= 0.2

dat <- cbind(colsumm, percentage, colsumm20)
dat <- as.data.frame(dat)

keep <- dat[which(dat[,"colsumm20"] == 1),]
keep_names <- rownames(keep)

  # harmonize between subset and keep
library(dplyr)
transposed_sub <- t(subset)
res <- subset(transposed_sub, rownames(transposed_sub) %in% keep_names)
res <- t(res)

# Put back the rest of the phenotypes
chunk1 <- phe[,1:35]

num <- ncol(res)
chunk2 <- res[,4:num]

merged <- cbind(chunk1, chunk2)
# Write out a CMM file and then use it for the normalize_phe.R code

#write.csv(merged, "pheno_cmm.csv", row.names = T)

