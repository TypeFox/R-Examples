### R code from vignette source 'PharmacoGx.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: setup (eval = FALSE)
###################################################
## options(keep.source=TRUE)


###################################################
### code chunk number 2: install-pkg (eval = FALSE)
###################################################
## source('http://bioconductor.org/biocLite.R')
## biocLite('PharmacoGx')


###################################################
### code chunk number 3: loadlib (eval = FALSE)
###################################################
## library(PharmacoGx)


###################################################
### code chunk number 4: download-example (eval = FALSE)
###################################################
## ## Example
## availablePSets()
## GDSC <- downloadPSet("GDSC") 


###################################################
### code chunk number 5: download-sig (eval = FALSE)
###################################################
## ## Example
## CMAP.sigs <- downloadPertSig("CMAP")


###################################################
### code chunk number 6: inconsistencies (eval = FALSE)
###################################################
##   library(Biobase)
##   library(PharmacoGx)
##   data("GDSCsmall")
##   data("CCLEsmall")
##   commonGenes <- intersect(featureNames(GDSCsmall, "rna"),
##                            featureNames(CCLEsmall,"rna"))
##   common <- intersectPSet(list('CCLE'=CCLEsmall,
##                                'GDSC'=GDSCsmall),
##                           intersectOn=c("cell.lines", "drugs"), strictIntersect=TRUE)
## 
##   
##   GDSC.auc <- summarizeSensitivityProfiles(
##                 common$GDSC,
##                 sensitivity.measure='auc_published', 
##                 summary.stat="median")
##   CCLE.auc <- summarizeSensitivityProfiles(
##                 common$CCLE,
##                 sensitivity.measure='auc_published', 
##                 summary.stat="median")
##   
##   GDSC.ic50 <- summarizeSensitivityProfiles(
##                 common$GDSC, 
##                 sensitivity.measure='ic50_published', 
##                 summary.stat="median")
##   CCLE.ic50 <- summarizeSensitivityProfiles(
##                 common$CCLE, 
##                 sensitivity.measure='ic50_published', 
##                 summary.stat="median")
##   
##   GDSCexpression <- summarizeMolecularProfiles(common$GDSC, 
##                                         cellNames(common$GDSC),
##                                         mDataType="rna",
##                                         features=commonGenes,
##                                         verbose=FALSE)
##   CCLEexpression <- summarizeMolecularProfiles(common$CCLE, 
##                                          cellNames(common$CCLE),
##                                          mDataType="rna",
##                                          features=commonGenes,
##                                          verbose=FALSE)
##   gg <- featureNames(common[[1]], 'rna')
##   cc <- cellNames(common[[1]])
##   
##   ge.cor <- sapply(cc, function (x, d1, d2) {
##     return (stats::cor(d1[ , x], d2[ , x], method="spearman",
##                 use="pairwise.complete.obs"))
##   }, d1=exprs(GDSCexpression), d2=exprs(CCLEexpression))
##   ic50.cor <- sapply(cc, function (x, d1, d2) {
##     return (stats::cor(d1[, x], d2[ , x], method="spearman",
##                 use="pairwise.complete.obs"))
##   }, d1=GDSC.ic50, d2=CCLE.ic50)
##   auc.cor <- sapply(cc, function (x, d1, d2) {
##     return (stats::cor(d1[ , x], d2[ , x], method="spearman",
##                 use="pairwise.complete.obs"))
##   }, d1=GDSC.auc, d2=CCLE.auc)
##   
##   w1 <- stats::wilcox.test(x=ge.cor, y=auc.cor,
##                            conf.int=TRUE, exact=FALSE)
##   w2 <- stats::wilcox.test(x=ge.cor, y=ic50.cor,
##                            conf.int=TRUE, exact=FALSE)
##   yylim <- c(-1, 1)
##   ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E",
##                 w1$p.value, w2$p.value)
##   boxplot(list("GE"=ge.cor,
##                "AUC"=auc.cor,
##                "IC50"=ic50.cor),
##           main="Concordance between cell lines",
##           ylab=expression(R[s]),
##           sub=ss,
##           ylim=yylim,
##           col="lightgrey",
##           pch=20,
##           border="black")
## 


###################################################
### code chunk number 7: fig2
###################################################

  library(Biobase)
  library(PharmacoGx)
  data("GDSCsmall")
  data("CCLEsmall")
  commonGenes <- intersect(featureNames(GDSCsmall, "rna"),
                           featureNames(CCLEsmall,"rna"))
  common <- intersectPSet(list('CCLE'=CCLEsmall,
                               'GDSC'=GDSCsmall),
                          intersectOn=c("cell.lines", "drugs"), strictIntersect=TRUE)

  
  GDSC.auc <- summarizeSensitivityProfiles(
                common$GDSC,
                sensitivity.measure='auc_published', 
                summary.stat="median")
  CCLE.auc <- summarizeSensitivityProfiles(
                common$CCLE,
                sensitivity.measure='auc_published', 
                summary.stat="median")
  
  GDSC.ic50 <- summarizeSensitivityProfiles(
                common$GDSC, 
                sensitivity.measure='ic50_published', 
                summary.stat="median")
  CCLE.ic50 <- summarizeSensitivityProfiles(
                common$CCLE, 
                sensitivity.measure='ic50_published', 
                summary.stat="median")
  
  GDSCexpression <- summarizeMolecularProfiles(common$GDSC, 
                                        cellNames(common$GDSC),
                                        mDataType="rna",
                                        features=commonGenes,
                                        verbose=FALSE)
  CCLEexpression <- summarizeMolecularProfiles(common$CCLE, 
                                         cellNames(common$CCLE),
                                         mDataType="rna",
                                         features=commonGenes,
                                         verbose=FALSE)
  gg <- featureNames(common[[1]], 'rna')
  cc <- cellNames(common[[1]])
  
  ge.cor <- sapply(cc, function (x, d1, d2) {
    return (stats::cor(d1[ , x], d2[ , x], method="spearman",
                use="pairwise.complete.obs"))
  }, d1=exprs(GDSCexpression), d2=exprs(CCLEexpression))
  ic50.cor <- sapply(cc, function (x, d1, d2) {
    return (stats::cor(d1[, x], d2[ , x], method="spearman",
                use="pairwise.complete.obs"))
  }, d1=GDSC.ic50, d2=CCLE.ic50)
  auc.cor <- sapply(cc, function (x, d1, d2) {
    return (stats::cor(d1[ , x], d2[ , x], method="spearman",
                use="pairwise.complete.obs"))
  }, d1=GDSC.auc, d2=CCLE.auc)
  
  w1 <- stats::wilcox.test(x=ge.cor, y=auc.cor,
                           conf.int=TRUE, exact=FALSE)
  w2 <- stats::wilcox.test(x=ge.cor, y=ic50.cor,
                           conf.int=TRUE, exact=FALSE)
  yylim <- c(-1, 1)
  ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E",
                w1$p.value, w2$p.value)
  boxplot(list("GE"=ge.cor,
               "AUC"=auc.cor,
               "IC50"=ic50.cor),
          main="Concordance between cell lines",
          ylab=expression(R[s]),
          sub=ss,
          ylim=yylim,
          col="lightgrey",
          pch=20,
          border="black")

    



###################################################
### code chunk number 8: load-cmap
###################################################
  library(PharmacoGx)
  require(xtable)
  data(CMAPsmall)
  drug.perturbation <- drugPerturbationSig(CMAPsmall, mDataType="rna")
  data(HDAC_genes)
  
  res <- apply(drug.perturbation[,,c("tstat", "fdr")],
               2, function(x, HDAC){ 
	    return(connectivityScore(x=x, 
	                             y=HDAC[,2,drop=FALSE], 
	                             method="gsea", nperm=100))
	}, HDAC=HDAC_genes)
  
  rownames(res) <- c("Connectivity", "P Value")
  res <- t(res)
  res <- res[order(res[,1], decreasing=TRUE),]
  xtable(res, 
    caption='Connectivity Score results for HDAC inhibitor gene signature.')


###################################################
### code chunk number 9: biomarkers
###################################################

  data(CCLEsmall)
  features <- featureNames(CCLEsmall, "rna")[
                          which(featureInfo(CCLEsmall,
                                            "rna")$Symbol == "NQO1")]
  sig.rna <- drugSensitivitySig(pSet=CCLEsmall, 
                            mDataType="rna", 
                            drugs=c("17-AAG"), 
                            features=features, 
                            sensitivity.measure="auc_published", 
                            molecular.summary.stat="median", 
                            sensitivity.summary.stat="median")
  sig.mut <- drugSensitivitySig(pSet=CCLEsmall, 
                            mDataType="mutation", 
                            drugs=c("PD-0325901"), 
                            features="BRAF", 
                            sensitivity.measure="auc_published", 
                            molecular.summary.stat="and", 
                            sensitivity.summary.stat="median")
  sig <- rbind(sig.rna, sig.mut)
  rownames(sig) <- c("17-AAG + NQO1","PD-0325901 + BRAF")
  colnames(sig) <- dimnames(sig.mut)[[3]]
  xtable(sig, caption='P Value of Gene-Drug Association')


###################################################
### code chunk number 10: sessionInfo
###################################################
toLatex(sessionInfo())


