### R code from vignette source 'cape.Rnw'

###################################################
### code chunk number 1: cape.Rnw:98-100
###################################################
library(cape)
data(obesity.cross)


###################################################
### code chunk number 2: cape.Rnw:135-136
###################################################
str(obesity.cross)


###################################################
### code chunk number 3: cape.Rnw:197-199
###################################################
obesity.cross <- select.pheno(obesity.cross, 
phenotypes = c("body_weight", "glucose", "insulin", "mom"))


###################################################
### code chunk number 4: cape.Rnw:207-208
###################################################
obesity.cross <- pheno2covar(obesity.cross, "mom")


###################################################
### code chunk number 5: cape.Rnw:236-237
###################################################
histPheno(obesity.cross)


###################################################
### code chunk number 6: cape.Rnw:250-251
###################################################
qqPheno(obesity.cross)


###################################################
### code chunk number 7: cape.Rnw:257-258
###################################################
plotPheno(obesity.cross, color.by = "mom", group.labels = c("non-obese", "obese"))


###################################################
### code chunk number 8: cape.Rnw:271-272
###################################################
obesity.cross <- norm.pheno(obesity.cross, mean.center = TRUE)


###################################################
### code chunk number 9: cape.Rnw:279-280
###################################################
histPheno(obesity.cross)


###################################################
### code chunk number 10: cape.Rnw:290-291
###################################################
qqPheno(obesity.cross)


###################################################
### code chunk number 11: cape.Rnw:324-325
###################################################
plotPhenoCor(obesity.cross, color.by = "mom", group.labels = c("non-obese", "obese"))


###################################################
### code chunk number 12: cape.Rnw:360-362
###################################################
obesity.cross <- get.eigentraits(obesity.cross, scale.pheno = FALSE,
normalize.pheno = FALSE)


###################################################
### code chunk number 13: cape.Rnw:374-377
###################################################
pdf("svd.pdf")
plotSVD(obesity.cross, orientation = "vertical")
dev.off()


###################################################
### code chunk number 14: cape.Rnw:417-418
###################################################
obesity.cross <- select.eigentraits(obesity.cross, traits.which = c(1,2))


###################################################
### code chunk number 15: cape.Rnw:476-480
###################################################
obesity.singlescan <- singlescan(obesity.cross, n.perm = 100, 
covar = "mom", scan.what = "eigentraits", alpha = c(0.01, 0.05), 
verbose = FALSE, use.kinship = FALSE, overwrite.alert = FALSE,
run.parallel = FALSE, n.cores = 2)


###################################################
### code chunk number 16: cape.Rnw:499-501
###################################################
plotSinglescan(data.obj = obesity.cross, singlescan.obj = obesity.singlescan, 
mark.chr = TRUE, mark.covar = FALSE)


###################################################
### code chunk number 17: cape.Rnw:552-554
###################################################
obesity.cross <- select.markers.for.pairscan(data.obj = obesity.cross, 
singlescan.obj = obesity.singlescan)


###################################################
### code chunk number 18: cape.Rnw:608-610
###################################################
plotSinglescan(data.obj = obesity.cross, singlescan.obj = obesity.singlescan, 
mark.chr = TRUE, show.rejected.markers = TRUE, standardized = TRUE)


###################################################
### code chunk number 19: cape.Rnw:618-621
###################################################
obesity.pairscan <- pairscan(data.obj = obesity.cross, covar = "mom", 
scan.what = "eigentraits", total.perm = 1000, min.per.genotype = 6, 
verbose = FALSE, overwrite.alert = FALSE, n.cores = 2)


###################################################
### code chunk number 20: cape.Rnw:687-689
###################################################
plotPairscan(data.obj = obesity.cross, pairscan.obj = obesity.pairscan, 
phenotype = "ET1", pdf.label = "Pairscan_Regression.pdf")


###################################################
### code chunk number 21: cape.Rnw:826-829
###################################################
obesity.cross <- error.prop(data.obj = obesity.cross, 
pairscan.obj = obesity.pairscan, 
perm = FALSE, verbose = FALSE, n.cores = 2)


###################################################
### code chunk number 22: cape.Rnw:840-843
###################################################
obesity.cross <- error.prop(data.obj = obesity.cross, 
pairscan.obj = obesity.pairscan, 
perm = TRUE, verbose = FALSE, n.cores = 2)


###################################################
### code chunk number 23: cape.Rnw:855-858
###################################################
obesity.cross <- calc.p(data.obj = obesity.cross, 
pairscan.obj = obesity.pairscan, pval.correction = "fdr",
n.cores = 2)


###################################################
### code chunk number 24: cape.Rnw:901-903
###################################################
obesity.cross <- direct.influence(data.obj = obesity.cross, 
pairscan.obj = obesity.pairscan, pval.correction = "fdr")


###################################################
### code chunk number 25: cape.Rnw:937-942
###################################################
pdf("variant_influences.pdf")
plotVariantInfluences(obesity.cross, p.or.q = 0.05, 
standardize = TRUE, not.tested.col = "lightgray", 
pheno.width = 8)
dev.off()


###################################################
### code chunk number 26: cape.Rnw:1027-1032
###################################################
pdf("Network_Collapsed.pdf")
obesity.cross <- get.network(obesity.cross, p.or.q = 0.05, 
collapse.linked.markers = TRUE, standardize = FALSE)
plotNetwork(obesity.cross, collapsed.net = TRUE)
dev.off()


###################################################
### code chunk number 27: cape.Rnw:1065-1070
###################################################
pdf("Network_Full.pdf")
obesity.cross <- get.network(obesity.cross, p.or.q = 0.05, 
collapse.linked.markers = FALSE, standardize = TRUE)
plotNetwork(obesity.cross, collapsed.net = FALSE)
dev.off()


###################################################
### code chunk number 28: cape.Rnw:1084-1087
###################################################
saveRDS(obesity.cross, "cross.RData")
saveRDS(obesity.singlescan, "crossSinglescan.RData")
saveRDS(obesity.pairscan, "crossPairscan.RData")


