# Linkage Analysis
# Cleaned for upload to Github 4/24/23 by ACS

# This is adapted from Arends Repository for Salvador, IJO, 2021

rm(list = ls())
source("/Volumes/Transcend\ 1/qtl_paper2_R/code/TexasMice2/functions.R")
setwd("/Volumes/Transcend\ 1/qtl_paper2_R")

## REALIZE THAT THE ORDER OF THE PHENOTYPE AND GENOTYPE FILES MATTERS TO MAPQTL FUNCTION
phenotypes <- read.csv("/Volumes/Transcend\ 1/qtl_paper2_R/data/raw/pathway_phe22_final_FINAL_final3.csv", row.names = 1)
genotypes <- read.table("data/gts.founder.txt", sep = "\t", row.names=1, colClasses="character")
map <- read.table("data/map.txt", sep = "\t", row.names=1)


# Scan starting CMM ASV data
newphenames <- colnames(phenotypes[181:ncol(phenotypes)]) # should start with fbfac7ca3133511d2bc1084a09656368
                                                          # should end with cbab5f434bd6c0322f7506d3416442eb

# Map additive and domiance deviation effects
#for(phe in newphenames){
#  res <- mapADQTL(numgeno, phenotypes, phe)
#  write.table(res, paste0("results/qtl_AD_", phe, ".txt"), sep = "\t", quote=FALSE)
#  png(paste0("results/qtl_AD_", phe, ".png"), width = 4 * 1024, height = 4 * 600, res = 300)
#    plotQTL(map, res, phe)
#  dev.off()
#}

# Map using marker as a factor (sum add and domdev)
for(phe in newphenames){
  res <- mapQTL(genotypes, phenotypes, phe)
  write.table(res, paste0("results/qtl_cmmASV_", phe, ".txt"), sep = "\t", quote=FALSE)
  #png(paste0("results/qtl_marker_", phe, ".png"), width = 4 * 1024, height = 4 * 600, res = 300)
    #plotQTL(map, res, phe, "marker")
  #dev.off()
}

## Interactive scans for Diet x QTL and Sex x QTL

for(phe in newphenames){
  res <- mapQTL_M_by_Diet(genotypes, phenotypes, phe)
  write.table(res, paste0("results/QxD_cmmASV_", phe, ".txt"), sep = "\t", quote=FALSE)
  #png(paste0("results/qtl_M_by_Diet_", phe, ".png"), width = 4 * 1024, height = 4 * 600, res = 300)
  #  plotQTL(map, res, phe, "marker")
  #dev.off()
}

for(phe in newphenames){
  res <- mapQTL_M_by_Sex(genotypes, phenotypes, phe)
  write.table(res, paste0("results/QxS_cmmASV_", phe, ".txt"), sep = "\t", quote=FALSE)
  #png(paste0("results/qtl_M_by_Sex_", phe, ".png"), width = 4 * 1024, height = 4 * 600, res = 300)
  #  plotQTL(map, res, phe, "marker")
  #dev.off()
}


