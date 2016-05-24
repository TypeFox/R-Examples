### R code from vignette source 'marray-notes.Rnw'

###################################################
### code chunk number 1: path-to-data
###################################################
library(DAAGbio)
path2data <- system.file("doc", package="DAAGbio")


###################################################
### code chunk number 2: load-limma
###################################################
library(limma)


###################################################
### code chunk number 3: readTargets
###################################################
targets <- readTargets("coralTargets.txt", path=path2data)
targets$FileName     # Display the file names


###################################################
### code chunk number 4: see-targets (eval = FALSE)
###################################################
## targets


###################################################
### code chunk number 5: read-images
###################################################
coralRG <- read.maimages(targets$FileName, source = "spot",
    path=path2data,
    other.columns=list(area="area", badspot="badspot"))


###################################################
### code chunk number 6: marray-notes.Rnw:205-206
###################################################
summary(coralRG$other$area)


###################################################
### code chunk number 7: marray-notes.Rnw:212-213
###################################################
plot(density(coralRG$other$area[,1]))


###################################################
### code chunk number 8: marray-notes.Rnw:225-228
###################################################
coralRG$genes <- readGAL(path=path2data)
coralRG$printer <- getLayout(coralRG$genes)
coralRG$printer


###################################################
### code chunk number 9: printseq
###################################################
plotprintseq()


###################################################
### code chunk number 10: marray-notes.Rnw:270-272
###################################################
spottypes<-readSpotTypes(path=path2data)
coralRG$genes$Status <- controlStatus(spottypes, coralRG)


###################################################
### code chunk number 11: marray-notes.Rnw:296-298
###################################################
imageplot(log2(coralRG$Rb[, 1]+1), layout = coralRG$printer,
          low="white", high="red")


###################################################
### code chunk number 12: marray-notes.Rnw:303-305
###################################################
imageplot(log2(coralRG$Rb[, 2]+1), layout = coralRG$printer,
          low="white", high="red")


###################################################
### code chunk number 13: marray-notes.Rnw:350-351
###################################################
plotMA(coralRG, array=1)


###################################################
### code chunk number 14: marray-notes.Rnw:368-370
###################################################
rawMA <- normalizeWithinArrays(coralRG, method = "none")
plotPrintTipLoess(rawMA, array=1)


###################################################
### code chunk number 15: marray-notes.Rnw:377-379
###################################################
MA <- normalizeWithinArrays(coralRG, method = "printtiploess")
plotPrintTipLoess(MA)


###################################################
### code chunk number 16: marray-notes.Rnw:382-383 (eval = FALSE)
###################################################
## boxplot(MA$M ~ col(MA$M), names = colnames(MA$M))


###################################################
### code chunk number 17: marray-notes.Rnw:390-392
###################################################
nMA <- normalizeBetweenArrays(MA)
boxplot(nMA$M ~ col(nMA$M), names = colnames(nMA$M))


###################################################
### code chunk number 18: marray-notes.Rnw:404-407
###################################################
wanted <- coralRG$genes$Status == "diff-exp ctl"
rawdeM <- rawMA$M[wanted, ]
pairs(rawdeM)


###################################################
### code chunk number 19: marray-notes.Rnw:415-418
###################################################
wanted <- coralRG$genes$Status == "diff-exp ctl"
deM <- nMA$M[wanted, ]
pairs(rawdeM)


###################################################
### code chunk number 20: marray-notes.Rnw:428-429
###################################################
imageplot(nMA$M[,5], layout=coralRG$printer)


###################################################
### code chunk number 21: marray-notes.Rnw:457-458
###################################################
design <- c(-1, 1, -1, 1, 0, 1)


###################################################
### code chunk number 22: marray-notes.Rnw:463-464
###################################################
fit <- lmFit(nMA, design)


###################################################
### code chunk number 23: marray-notes.Rnw:496-499
###################################################
efit <- eBayes(fit)
qqt(efit$t, df = efit$df.prior + efit$df.residual, pch = 16,
    cex = 0.2)


###################################################
### code chunk number 24: marray-notes.Rnw:511-514
###################################################
options(digits = 3)
topvals <- topTable(efit, number = 50)
topvals


###################################################
### code chunk number 25: marray-notes.Rnw:526-532
###################################################
plot(efit$coef, efit$lods, pch = 16, cex = 0.2, xlab = "log(fold change)",
    ylab = "log(odds)")
ord <- order(efit$lods, decreasing = TRUE)
top8 <- ord[1:8]
text(efit$coef[top8], efit$lods[top8], labels = coralRG$genes[top8,
    "Name"], cex = 0.8, col = "blue")


###################################################
### code chunk number 26: marray-notes.Rnw:556-562
###################################################
coral2RG <- read.maimages(targets$FileName,
                         source = "spot",
                         path=path2data,
                         wt.fun=wtarea(100))
coral2RG$genes <- readGAL(path=path2data)
coral2RG$printer <- getLayout(coral2RG$genes)


###################################################
### code chunk number 27: marray-notes.Rnw:569-571
###################################################
MA2 <- normalizeWithinArrays(coral2RG, method = "printtiploess")
plotPrintTipLoess(MA2)


###################################################
### code chunk number 28: marray-notes.Rnw:574-577
###################################################
boxplot(MA2$M ~ col(MA2$M), names = colnames(MA2$M))
nMA2 <- normalizeBetweenArrays(MA2)
boxplot(nMA2$M ~ col(nMA2$M), names = colnames(nMA2$M))


###################################################
### code chunk number 29: marray-notes.Rnw:583-584
###################################################
imageplot(nMA2$M[,5], layout=coral2RG$printer)


###################################################
### code chunk number 30: marray-notes.Rnw:590-592
###################################################
design <- c(-1, 1, -1, 1, 0, 1)
fit2 <- lmFit(nMA2, design)


###################################################
### code chunk number 31: marray-notes.Rnw:597-600
###################################################
efit2 <- eBayes(fit2)
qqt(efit2$t, df = efit2$df.prior + efit2$df.residual, pch = 16,
    cex = 0.2)


###################################################
### code chunk number 32: marray-notes.Rnw:605-607
###################################################
options(digits = 3)
topTable(efit2, number = 50)


###################################################
### code chunk number 33: marray-notes.Rnw:612-615
###################################################
## Get & store results with & without weights
topvals2 <- topTable(efit2, number = 50)
cbind(row.names(topvals), row.names(topvals2))


###################################################
### code chunk number 34: marray-notes.Rnw:618-619
###################################################
sum(row.names(topvals)%in%row.names(topvals2))


###################################################
### code chunk number 35: marray-notes.Rnw:646-647
###################################################
imgplot(coralRG$R[, 1], layout = coralRG$printer)


