# Defined Functions
# Cleaned for GitHub 4/24/23 by ACS

# Adapted from Arends repository for Salvador, IJO, 2021

chrs <- c(1:19, "X", "Y")

####################################################################################################
# Interactive QTL 
mapQTL_M_by_Diet <- function(genotypes, phenotypes, phe = "fat_final"){
  pvals <- c()
  for(x in 1:nrow(genotypes)){
    pheno <- phenotypes[, phe]
    sex <- phenotypes[, "sex"]
    diet <- phenotypes[, "diet"]
    marker <- as.factor(genotypes[x,])
    mymodel <- lm(pheno ~ sex + diet + marker + marker:diet)
    pval.sex <- anova(mymodel)["sex", "Pr(>F)"]
    pval.diet <- anova(mymodel)["diet", "Pr(>F)"]
    pval <- anova(mymodel)["marker", "Pr(>F)"]
    pval.dietmarker <- anova(mymodel)["diet:marker", "Pr(>F)"]
    
    
    dH <- as.numeric(coefficients(mymodel)["dietWestern:markerH"])
    dF <- as.numeric(coefficients(mymodel)["dietWestern:markerF"])
    varExplained <- round(100 * anova(mymodel)["diet:marker", "Sum Sq"] / sum(anova(mymodel)[, "Sum Sq"]),2)
    
    pvals <- rbind(pvals, c(pval.sex, pval.diet, pval, pval.dietmarker, dH, dF, varExplained))
  }
  colnames(pvals) <- c("sex", "diet", "marker", "diet:marker", "dH", "dF", "varExplained")
  rownames(pvals) <- rownames(genotypes)
  return(pvals)
}

mapQTL_M_by_Sex <- function(genotypes, phenotypes, phe = "fat_final"){
  pvals <- c()
  for(x in 1:nrow(genotypes)){
    pheno <- phenotypes[, phe]
    sex <- phenotypes[, "sex"]
    diet <- phenotypes[, "diet"]
    marker <- as.factor(genotypes[x,])
    mymodel <- lm(pheno ~ diet + sex + marker + marker:sex) #Version 7
#    mymodel <- lm(pheno ~ sex + marker + marker:sex) # Version 6
    pval.sex <- anova(mymodel)["sex", "Pr(>F)"]
    pval.diet <- anova(mymodel)["diet", "Pr(>F)"]
    pval <- anova(mymodel)["marker", "Pr(>F)"]
    pval.sexmarker <- anova(mymodel)["sex:marker", "Pr(>F)"]
    
    
    dH <- as.numeric(coefficients(mymodel)["sexM:markerH"])
    dF <- as.numeric(coefficients(mymodel)["sexM:markerF"])
    varExplained <- round(100 * anova(mymodel)["sex:marker", "Sum Sq"] / sum(anova(mymodel)[, "Sum Sq"]),2)
    
#    pvals <- rbind(pvals, c(pval.sex, pval.diet, pval, pval.sexmarker, dH, dF, varExplained))
    pvals <- rbind(pvals, c(pval.sex, pval.diet, pval, pval.sexmarker, dH, dF, varExplained))
  }
  colnames(pvals) <- c("sex", "diet", "marker", "sex:marker", "dH", "dF", "varExplained")
  rownames(pvals) <- rownames(genotypes)
  return(pvals)
}

mapQTL_M_by_Sex_DietSpecific <- function(genotypes, phenotypes, phe = "fat_final"){
  pvals <- c()
  for(x in 1:nrow(genotypes)){
    pheno <- phenotypes[, phe]
    sex <- phenotypes[, "sex"]
    #    diet <- phenotypes[, "diet"]
    marker <- as.factor(genotypes[x,])
    #    mymodel <- lm(pheno ~ sex + diet + marker + marker:sex)
    mymodel <- lm(pheno ~ sex + marker + marker:sex)
    pval.sex <- anova(mymodel)["sex", "Pr(>F)"]
    #    pval.diet <- anova(mymodel)["diet", "Pr(>F)"]
    pval <- anova(mymodel)["marker", "Pr(>F)"]
    pval.sexmarker <- anova(mymodel)["sex:marker", "Pr(>F)"]
    
    
    dH <- as.numeric(coefficients(mymodel)["markerH"])
    dF <- as.numeric(coefficients(mymodel)["markerF"])
    varExplained <- round(100 * anova(mymodel)["marker", "Sum Sq"] / sum(anova(mymodel)[, "Sum Sq"]),2)
    
    #    pvals <- rbind(pvals, c(pval.sex, pval.diet, pval, pval.sexmarker, dH, dF, varExplained))
    pvals <- rbind(pvals, c(pval.sex, pval, pval.sexmarker, dH, dF, varExplained))
  }
  colnames(pvals) <- c("sex", "marker", "sex:marker", "dH", "dF", "varExplained")
  rownames(pvals) <- rownames(genotypes)
  return(pvals)
}
###########################################################################################################

getPositions <- function(map, gap = 25000000){
  chr.start <- 0
  nmap <- c()
  for(chr in chrs){
    positions <- chr.start + map[which(map[, "chr"] == chr), "bp_mm10"]
    chr.start <- gap + max(positions, na.rm = TRUE)
    nmap <- rbind(nmap, cbind(chr = rep(chr, length(positions)), pos = positions))
  }
  rownames(nmap) <- rownames(map)
  return(nmap)
}

getChrSE <- function(map, gap = 25000000){
  chr.start <- 0
  nmap <- c()
  for(chr in chrs){
    positions <- chr.start + map[which(map[, "chr"] == chr), "bp_mm10"]
    chr.end <- max(positions, na.rm = TRUE)
    nmap <- rbind(nmap, cbind(chr.start, chr.end))
    chr.start <- gap + max(positions, na.rm = TRUE)
  }
  rownames(nmap) <- chrs
  return(nmap)
}

getChrIdx <- function(map){
  chrstart <- 0
  nmap <- c()
  for(chr in chrs){
    nmar <- chrstart + length(which(map[, "chr"] == chr))
    nmap <- c(nmap, nmar)
    chrstart <- nmar
  }
  return(nmap)
}

plotQTL <- function(map, res, phe, what = c("add", "dom"), what2 = c("", ""), nTests = 500, colz = c("black", "blue", "magenta3", "orange", "forestgreen"), plotChr = NULL){ 
#plotQTL <- function(map, res, phe, what = c("add", "dom"), what2 = c("", ""), nTests = 500, colz = c("black", "forestgreen", "palevioletred4", "orange", "forestgreen"), plotChr = NULL){ #forLDL
  nmap <- getPositions(map, 50000000)
  cmap <- getChrSE(map, 50000000)
  maxY <- 0
  for(y in what){
    maxY <- max(maxY, 1.25 * max(c(-log10(abs(res[, y]))), na.rm = TRUE))
  }
  if(!is.null(plotChr)){
    submap <- nmap[which(nmap[,1] == plotChr),]
    plot(c(min(as.numeric(submap[, "pos"])),max(as.numeric(submap[, "pos"]))), c(0,maxY), 
         t = 'n', xlab="Position (Mb)", ylab = "-log10(p)", xaxt = "n", xaxs = 'i', yaxs = 'i', 
         las = 2, main = phe)  
    rect(105000000, -1, 135000000, 20, col = rgb(1,0.75,0.78,0.5), border = NA) #for HDL
    #rect(135000000, -1, 153000000, 20, col= rgb(1,0.75,0.78,0.5), border = NA ) #for LDL
    rect(153000000, -1, 183000000, 20, col = rgb(0.8,0.65,1.0,0.5), border = NA) #for HDL
    #rect(153000000, -1, 188000000, 20, col = rgb(0.8,0.65,1.0,0.5), border = NA)
  }else{
    #    plot(c(0,max(as.numeric(nmap[, "pos"]))), c(-maxY,maxY), t = 'n', cex.main = 1.5, cex.axis= 1.5, cex.lab = 1.5, xlab="Chr", ylab = "LOD", # If you want to show the direction of the effects
    plot(c(0,max(as.numeric(nmap[, "pos"]))), c(0,maxY), t = 'n', cex.main = 1.5, cex.axis= 1.5, cex.lab = 1.5, xlab="Chr", ylab = "LOD", 
         xaxt = "n", xaxs = 'i', las = 2, main = phe)
  }
  abline(v = cmap[,1] - 25000000, col = "lightgray")
  #abline(v = cmap[,2], col = "lightgray")
  cnt <- 1
  #colz <- c("black", "blue", "magenta3", "orange", "forestgreen")
  for(chr in chrs){
    onChr <- which(nmap[, "chr"] == chr)
    positions <- as.numeric(nmap[onChr, "pos"])
    ccnt <- 1
    for(y in what){
      LOD <- -log10(abs(res[onChr, y]))
#      LOD <- LOD * sign(res[onChr, y])
      points(positions, LOD, t = 'l', col=colz[ccnt], lwd=2)
      ccnt <- ccnt + 1
    }
    cnt <- cnt + 1
  }
  abline(h = -log10(0.01 / nTests), col = "green", lwd=2, lty=2)
  abline(h = -log10(0.05 / nTests), col = "orange", lwd=2, lty=2)
  #abline(h = -log10(0.1 / nTests), col = "red", lwd=2, lty=2)
  abline(h = 0, col = "black")
  #abline(h = log10(0.1 / nTests), col = "red", lwd=2, lty=2)
  abline(h = log10(0.05 / nTests), col = "orange", lwd=2, lty=2)
  abline(h = log10(0.01 / nTests), col = "green", lwd=2, lty=2)
  if(!is.null(plotChr)){
    submap <- nmap[which(nmap[,1] == plotChr),]
    axis(1, at = seq(0,max(as.numeric(submap[, "pos"])), 10000000), seq(0,max(as.numeric(submap[, "pos"])), 10000000) / 1000000)
  }else{
    axis(1, at = cmap[, "chr.start"] + ((cmap[, "chr.end"] - cmap[, "chr.start"]) / 2.0), rownames(getChrSE(map)))
  }
  cat(what2, "\n")
  legend("topright", what2, col=colz[1:length(what)], lwd=2, bg="white")
  #legend("topleft", c("5% FDR", "1% FDR"), col=c("orange","green"), lwd=2, lty=2, bg="white")
  #legend("topleft", c("10% FDR", "5% FDR", "1% FDR"), col=c("red", "orange","green"), lwd=2, lty=2, bg="white")
}

mapADQTL <- function(genotypes, phenotypes, phe = "fat_final"){
  pvals <- c()
  for(x in 1:nrow(genotypes)){
    pheno <- phenotypes[, phe]
    sex <- phenotypes[, "sex"]
    diet <- phenotypes[, "diet"]
    marker.add <- genotypes[x,]
    marker.dom <- as.numeric(genotypes[x,] == 0)
    mymodel <- lm(pheno ~ sex + diet + marker.add + marker.dom)
    pval.sex <- anova(mymodel)["sex", "Pr(>F)"]
    pval.diet <- anova(mymodel)["diet", "Pr(>F)"]
    pval.add <- anova(mymodel)["marker.add", "Pr(>F)"]
    pval.dom <- anova(mymodel)["marker.dom", "Pr(>F)"]
    pvals <- rbind(pvals, c(pval.sex, pval.diet, pval.add, pval.dom))
  }
  colnames(pvals) <- c("sex", "diet", "add", "dom")
  rownames(pvals) <- rownames(genotypes)
  return(pvals)
}

mapQTL <- function(genotypes, phenotypes, phe = "fat_final"){
  pvals <- c()
  for(x in 1:nrow(genotypes)){
    pheno <- phenotypes[, phe] #REALIZE THAT THE ORDER OF THE PHENOTYPE AND GENOTYPE FILE MATTERS
    sex <- phenotypes[, "sex"]#SAMPLES ARE DE-IDENTIFIED DURING THESE STEPS!!!
    diet <- phenotypes[, "diet"]
    marker <- as.factor(genotypes[x,])
    mymodel <- lm(pheno ~ sex * diet + marker)
    pval.sex <- anova(mymodel)["sex", "Pr(>F)"]
    pval.diet <- anova(mymodel)["diet", "Pr(>F)"]
    pval.sexdiet <- anova(mymodel)["sex:diet", "Pr(>F)"]
    pval <- anova(mymodel)["marker", "Pr(>F)"]
    
    dH <- as.numeric(coefficients(mymodel)["markerH"])
    dF <- as.numeric(coefficients(mymodel)["markerF"])
    varExplained <- round(100 * anova(mymodel)["marker", "Sum Sq"] / sum(anova(mymodel)[, "Sum Sq"]),2)
    
    pvals <- rbind(pvals, c(pval.sex, pval.diet, pval.sexdiet, pval, dH, dF, varExplained))
  }
  colnames(pvals) <- c("sex", "diet", "sex:diet", "marker", "dH", "dF", "varExplained")
  rownames(pvals) <- rownames(genotypes)
  return(pvals)
}

mapQTLcim <- function(genotypes, phenotypes, phe = "fat_gain", mname = "backupUNC050383757"){
  pvals <- c()
  topmarker <- as.factor(genotypes[mname,])
  for(x in 1:nrow(genotypes)){
    pheno <- phenotypes[, phe]
    sex <- phenotypes[, "sex"]
    diet <- phenotypes[, "diet"]
    marker <- as.factor(genotypes[x,])
    
    mymodel <- lm(pheno ~ sex * diet + topmarker + marker)
    pval.sex <- anova(mymodel)["sex", "Pr(>F)"]
    pval.diet <- anova(mymodel)["diet", "Pr(>F)"]
    pval.sexdiet <- anova(mymodel)["sex:diet", "Pr(>F)"]
    pval.topmarker <- anova(mymodel)["topmarker", "Pr(>F)"]
    pval <- anova(mymodel)["marker", "Pr(>F)"]
    
    dH <- as.numeric(coefficients(mymodel)["markerH"])
    dF <- as.numeric(coefficients(mymodel)["markerF"])
    varExplained <- round(100 * anova(mymodel)["marker", "Sum Sq"] / sum(anova(mymodel)[, "Sum Sq"]),2)
    
    pvals <- rbind(pvals, c(pval.sex, pval.diet, pval.sexdiet, pval.topmarker, pval, dH, dF, varExplained))
  }
  colnames(pvals) <- c("sex", "diet", "sex:diet", mname, "marker", "dH", "dF", "varExplained")
  rownames(pvals) <- rownames(genotypes)
  return(pvals)
}

mapQTLdiet <- function(genotypes, phenotypes, diet = "Keto", phe = "fat_final"){
  ondiet <- which(phenotypes[, "diet"]==diet)
  phenosubset <- phenotypes[ondiet,]
  genosubset <- genotypes[,ondiet]
  
  pvals <- c()
  for(x in 1:nrow(genosubset)){
    pheno <- phenosubset[, phe]
    sex <- phenosubset[, "sex"]
    marker <- as.factor(genosubset[x,])
    mymodel <- lm(pheno ~ sex + marker)
    pval.sex <- anova(mymodel)["sex", "Pr(>F)"]
    pval <- anova(mymodel)["marker", "Pr(>F)"]
    
    dH <- as.numeric(coefficients(mymodel)["markerH"])
    dF <- as.numeric(coefficients(mymodel)["markerF"])
    varExplained <- round(100 * anova(mymodel)["marker", "Sum Sq"] / sum(anova(mymodel)[, "Sum Sq"]),2)

    pvals <- rbind(pvals, c(pval.sex, pval, dH, dF, varExplained))
  }
  colnames(pvals) <- c("sex", "marker", "dH", "dF", "varExplained")
  rownames(pvals) <- rownames(genosubset)
  return(pvals)
}

mapQTLsex <- function(genotypes, phenotypes, sex = "M", phe = "fat_final"){
  ondiet <- which(phenotypes[, "sex"]==sex)
  phenosubset <- phenotypes[ondiet,]
  genosubset <- genotypes[,ondiet]
  
  pvals <- c()
  for(x in 1:nrow(genosubset)){
    pheno <- phenosubset[, phe]
    diet <- phenosubset[, "diet"]
    marker <- as.factor(genosubset[x,])
    if(length(table(marker)) != 1){
      mymodel <- lm(pheno ~ diet + marker)
      pval.diet <- anova(mymodel)["diet", "Pr(>F)"]
      pval <- anova(mymodel)["marker", "Pr(>F)"]
      dH <- as.numeric(coefficients(mymodel)["markerH"])
      dF <- as.numeric(coefficients(mymodel)["markerF"])
      varExplained <- round(100 * anova(mymodel)["marker", "Sum Sq"] / sum(anova(mymodel)[, "Sum Sq"]),2)

      pvals <- rbind(pvals, c(pval.diet, pval, dH, dF, varExplained))
    }else{
      pvals <- rbind(pvals, c(NA, NA, NA, NA, NA))
    }
  }
  colnames(pvals) <- c("diet", "marker", "dH", "dF", "varExplained")
  rownames(pvals) <- rownames(genosubset)
  return(pvals)
}

mapQTL_DietBysex <- function(genotypes, phenotypes, sex = "M", phe = "fat_final"){
  ondiet <- which(phenotypes[, "sex"]==sex)
  phenosubset <- phenotypes[ondiet,]
  genosubset <- genotypes[,ondiet]
  
  pvals <- c()
  for(x in 1:nrow(genosubset)){
    pheno <- phenosubset[, phe]
    diet <- phenosubset[, "diet"]
    marker <- as.factor(genosubset[x,])
    if(length(table(marker)) != 1){
      mymodel <- lm(pheno ~ diet + marker + diet:marker)
      pval.diet <- anova(mymodel)["diet", "Pr(>F)"]
      pval.marker <- anova(mymodel)["marker", "Pr(>F)"]
      pval.int <- anova(mymodel)["diet:marker", "Pr(>F)"]
      dH <- as.numeric(coefficients(mymodel)["markerH"])
      dF <- as.numeric(coefficients(mymodel)["markerF"])
      varExplained <- round(100 * anova(mymodel)["marker", "Sum Sq"] / sum(anova(mymodel)[, "Sum Sq"]),2)

      pvals <- rbind(pvals, c(pval.diet, pval.marker, pval.int, dH, dF, varExplained))
    }else{
      pvals <- rbind(pvals, c(NA, NA, NA, NA, NA, NA))
    }
  }
  colnames(pvals) <- c("diet", "marker", "diet:marker", "dH", "dF", "varExplained")
  rownames(pvals) <- rownames(genosubset)
  return(pvals)
}

mapQTLsexdiet <- function(genotypes, phenotypes, sex = "M", diet = "Keto", phe = "fat_final"){
  ondiet <- which(phenotypes[, "sex"]==sex & phenotypes[, "diet"]==diet)
  phenosubset <- phenotypes[ondiet,]
  genosubset <- genotypes[,ondiet]
  
  pvals <- c()
  for(x in 1:nrow(genosubset)){
    pheno <- phenosubset[, phe]
    marker <- as.factor(genosubset[x,])
    if(length(table(marker)) != 1){
      mymodel <- lm(pheno ~ marker)
      pval <- anova(mymodel)["marker", "Pr(>F)"]
      dH <- as.numeric(coefficients(mymodel)["markerH"])
      dF <- as.numeric(coefficients(mymodel)["markerF"])
      varExplained <- round(100 * anova(mymodel)["marker", "Sum Sq"] / sum(anova(mymodel)[, "Sum Sq"]),2)

      pvals <- rbind(pvals, c(pval, dH, dF, varExplained))
    }else{
      pvals <- rbind(pvals, c(NA, NA, NA, NA))
    }
  }
  colnames(pvals) <- c("marker", "dH", "dF", "varExplained")
  rownames(pvals) <- rownames(genosubset)
  return(pvals)
}

getpeaks <- function(qtlprofiles, cutoff = 3.0){
  #cat("Starting peak detection above",cutoff, "nrow", nrow(qtlprofiles),"\n")
  peaks <- vector("list", nrow(qtlprofiles))
  for(x in 1:nrow(qtlprofiles)){
    peak <- FALSE
    curmax <- 0
    curmaxindex <- 1
    marker <- 1
    maximums <- NULL
    for(ab in (qtlprofiles[x,]>cutoff | qtlprofiles[x,]<(-cutoff))){
      if(ab){
        peak <- TRUE
        if(qtlprofiles[x,marker]/abs(qtlprofiles[x,marker]) > 0){
          if(qtlprofiles[x,marker] > curmax){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- marker
          }
        }else{
          if(qtlprofiles[x,marker] < (-curmax)){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- -marker
          }
        }
        if(ncol(qtlprofiles)==marker){
          if(curmax!=0) maximums <- c(maximums,curmaxindex)
        }
      }else{
        if(curmax!=0) maximums <- c(maximums,curmaxindex)
        peak <- FALSE
        curmax <- 0
      }
      marker <- marker+1
    }
    if(length(maximums) > 0) peaks[[x]] <- maximums
  }
  return(peaks)
}