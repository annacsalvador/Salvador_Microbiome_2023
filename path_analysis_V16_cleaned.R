# SEM Path Analysis
# Cleaned for GitHub 4/24/23 by ACS

rm(list=ls())

# Based on https://www.rpubs.com/tbihansk/302732
#install.packages("OpenMx")
#install.packages("lavaan")
#install.packages("semPlot")
#install.packages("lme4")
#install.packages("tidyverse")
#install.packages("xml2")
#install.packages("knitr")
#install.packages("kableExtra")
#install.packages("GGally")

library(lavaan)
library(semPlot)
library(OpenMx)
library(tidyverse)
library(knitr)
library(kableExtra)
library(GGally)

#### Pathway analysis!
setwd("/Volumes/Transcend\ 1/Active_Projects/qtl_paper2_R/")

# Requires genotypes and phenotypes to be in the same file
phenotypes <- read.csv("data/raw/pathway_phe22_final_FINAL.csv", row.names = 1)

#genotypes <- read.table("data/gts.founder.txt", sep = "\t", row.names=1, colClasses="character")

# ASV names
# ecfbbfd673f8a3e1c34718d968e6b069	# Uncultured Rikenella # Asvq7
# b89e2f23a606787577644904f2b6aed3 # Ruminiclostridium 9 # Asvq16
# X18f60eae2c3a26bab8efa0f9823e31b4 # Uncultured Bilophila # Asvq17 # remember this phe name has a leading X

# d6554036e1a2bd9f0586972111720879 # Uncultured Rikenelleceae RC9 Gut Group # Asvq23
# a371498f7ad6d1368940364ff67938b6 # Parabacteroides Goldsteini CL02T12C30 # Asvq24

# other phe names
# wunifrac_Axis.1 # wufpc1q
# bray_Axis.2 # bcpc2q_jpc2q
# jaccard_Axis.2 # bcpc2q_jpc2q
# fat_gain # fmgq1
# hdl_avg # hdlq1



# DWLS is only confirmed to be appropriate for ordinal

library(dplyr) # Collapse into ordinal quantiles, use histograms for sanity check. 

hist(phenotypes$X18f60eae2c3a26bab8efa0f9823e31b4)
phenotypes$X18f60eae2c3a26bab8efa0f9823e31b4 <- ntile(phenotypes$X18f60eae2c3a26bab8efa0f9823e31b4, 4)
hist(phenotypes$X18f60eae2c3a26bab8efa0f9823e31b4)

hist(phenotypes$b89e2f23a606787577644904f2b6aed3)
phenotypes$b89e2f23a606787577644904f2b6aed3 <- ntile(phenotypes$b89e2f23a606787577644904f2b6aed3, 4)
hist(phenotypes$b89e2f23a606787577644904f2b6aed3)

hist(phenotypes$ecfbbfd673f8a3e1c34718d968e6b069)
phenotypes$ecfbbfd673f8a3e1c34718d968e6b069 <- ntile(phenotypes$ecfbbfd673f8a3e1c34718d968e6b069, 4)
hist(phenotypes$ecfbbfd673f8a3e1c34718d968e6b069)

hist(phenotypes$d6554036e1a2bd9f0586972111720879)
phenotypes$d6554036e1a2bd9f0586972111720879 <- ntile(phenotypes$d6554036e1a2bd9f0586972111720879, 4)
hist(phenotypes$d6554036e1a2bd9f0586972111720879)

hist(phenotypes$a371498f7ad6d1368940364ff67938b6)
phenotypes$a371498f7ad6d1368940364ff67938b6 <- ntile(phenotypes$a371498f7ad6d1368940364ff67938b6, 4)
hist(phenotypes$a371498f7ad6d1368940364ff67938b6)

hist(phenotypes$bray_Axis.2)
phenotypes$bray_Axis.2 <- ntile(phenotypes$bray_Axis.2, 4)
hist(phenotypes$bray_Axis.2)

hist(phenotypes$wunifrac_Axis.1)
phenotypes$wunifrac_Axis.1 <- ntile(phenotypes$wunifrac_Axis.1, 4)
hist(phenotypes$wunifrac_Axis.1)


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
# Test the causal models Chr 1

sem1d <- 'bray_Axis.2 ~ diet + bcpc2q_jpc2q + ecfbbfd673f8a3e1c34718d968e6b069 + b89e2f23a606787577644904f2b6aed3 + X18f60eae2c3a26bab8efa0f9823e31b4
          X18f60eae2c3a26bab8efa0f9823e31b4 ~ bcpc2q_jpc2q
          b89e2f23a606787577644904f2b6aed3 ~ bcpc2q_jpc2q
          X18f60eae2c3a26bab8efa0f9823e31b4 ~~ b89e2f23a606787577644904f2b6aed3
          hdl_avg ~ sex + diet + bcpc2q_jpc2q
          fat_gain ~ sex + bcpc2q_jpc2q
' 

fit_data1d <- cfa(sem1d, data = phenotypes, 
                  estimator = "DWLS")
summary(fit_data1d, fit.measures = TRUE, standardized=T,rsquare=T)
semPaths(fit_data1d,"std",layout = 'tree', edge.label.cex=.9, curvePivot = TRUE)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
# Refine sem1d

sem1d.R <- 'bray_Axis.2 ~ diet + ecfbbfd673f8a3e1c34718d968e6b069 + b89e2f23a606787577644904f2b6aed3 + X18f60eae2c3a26bab8efa0f9823e31b4
            X18f60eae2c3a26bab8efa0f9823e31b4 ~ bcpc2q_jpc2q
            b89e2f23a606787577644904f2b6aed3 ~ bcpc2q_jpc2q
            X18f60eae2c3a26bab8efa0f9823e31b4 ~~ b89e2f23a606787577644904f2b6aed3
            hdl_avg ~ sex + diet + bcpc2q_jpc2q
            fat_gain ~ sex + bcpc2q_jpc2q
' 
fit_data1d.R <- cfa(sem1d.R, data = phenotypes, 
                    estimator = "DWLS")
summary(fit_data1d.R, fit.measures = TRUE, standardized=T,rsquare=T)
semPaths(fit_data1d.R,"std",layout = 'tree', edge.label.cex=.9, curvePivot = TRUE)
####################################################################################################

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#*# Test the causal models Chr 16

sem2c <- 'wunifrac_Axis.1 ~ diet + wufpc1q + wufpc1q*diet + d6554036e1a2bd9f0586972111720879 + a371498f7ad6d1368940364ff67938b6
          a371498f7ad6d1368940364ff67938b6 ~ wufpc1q*diet
          d6554036e1a2bd9f0586972111720879 ~ wufpc1q*diet
          a371498f7ad6d1368940364ff67938b6 ~~ d6554036e1a2bd9f0586972111720879
'

fit_data2c <- cfa(sem2c, data = phenotypes, 
                  estimator = "DWLS")
summary(fit_data2c, fit.measures = TRUE, standardized=T,rsquare=T)
semPaths(fit_data2c,"std",layout = 'tree', edge.label.cex=.9, curvePivot = TRUE)


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
# Refine for DWLS

sem2c.Rb <- 'wunifrac_Axis.1 ~  diet + a371498f7ad6d1368940364ff67938b6
              a371498f7ad6d1368940364ff67938b6 ~ wufpc1q*diet
              d6554036e1a2bd9f0586972111720879 ~ wufpc1q*diet
              a371498f7ad6d1368940364ff67938b6 ~~ d6554036e1a2bd9f0586972111720879
'

fit_data2c.Rb <- cfa(sem2c.Rb, data = phenotypes, 
                  estimator = "DWLS")
summary(fit_data2c.Rb, fit.measures = TRUE, standardized=T,rsquare=T)
semPaths(fit_data2c.Rb,"std",layout = 'tree', edge.label.cex=.9, curvePivot = TRUE)


