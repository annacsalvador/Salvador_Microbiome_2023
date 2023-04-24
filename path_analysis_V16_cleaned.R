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
setwd("/Volumes/Transcend1/Active_Projects/qtl_paper2_R/")

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

sem1e <- 'bray_Axis.2 ~ diet + bcpc2q_jpc2q
          X18f60eae2c3a26bab8efa0f9823e31b4 ~ diet + bcpc2q_jpc2q
          b89e2f23a606787577644904f2b6aed3 ~ sex + diet + sex*diet + bcpc2q_jpc2q
          ecfbbfd673f8a3e1c34718d968e6b069 ~ sex + diet + bcpc2q_jpc2q
          X18f60eae2c3a26bab8efa0f9823e31b4 ~~ b89e2f23a606787577644904f2b6aed3
          hdl_avg ~ sex + diet + bcpc2q_jpc2q
          fat_gain ~ sex + bcpc2q_jpc2q
' 

fit_data1e <- cfa(sem1e, data = phenotypes, 
                  estimator = "DWLS")
summary(fit_data1e, fit.measures = TRUE, standardized=T,rsquare=T)
semPaths(fit_data1e,"std",layout = 'tree', edge.label.cex=.9, curvePivot = TRUE)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
# Refine sem1e

sem1e.R <- 'bray_Axis.2 ~ diet + bcpc2q_jpc2q
          X18f60eae2c3a26bab8efa0f9823e31b4 ~ diet + bcpc2q_jpc2q
          b89e2f23a606787577644904f2b6aed3 ~ diet
          ecfbbfd673f8a3e1c34718d968e6b069 ~ bcpc2q_jpc2q
          X18f60eae2c3a26bab8efa0f9823e31b4 ~~ b89e2f23a606787577644904f2b6aed3
          hdl_avg ~ sex + bcpc2q_jpc2q
          fat_gain ~ sex + bcpc2q_jpc2q
' 

fit_data1e.R <- cfa(sem1e.R, data = phenotypes, 
                  estimator = "DWLS")
summary(fit_data1e.R, fit.measures = TRUE, standardized=T,rsquare=T)
semPaths(fit_data1e.R,"std",layout = 'tree', edge.label.cex=.9, curvePivot = TRUE)
      # Note that the most assertive arguments we would make about sem1d are detected here as covariate paths.
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


