### R code from vignette source 'RFLPtools.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: load
###################################################
library(RFLPtools)


###################################################
### code chunk number 2: RFLPqc
###################################################
Dir <- system.file("extdata", package = "RFLPtools") # input directory 
filename <- file.path(Dir, "AZ091016_report.txt")
RFLP1 <- read.rflp(file = filename)
str(RFLP1)

RFLP2 <- RFLPqc(RFLP1, rm.band1 = FALSE) # identical to RFLP1
identical(RFLP1, RFLP2)

RFLP3 <- RFLPqc(RFLP1)
str(RFLP3)

RFLP4 <- RFLPqc(RFLP1, rm.band1 = TRUE, QC.rm = TRUE)
str(RFLP4)


###################################################
### code chunk number 3: eucl
###################################################
data(RFLPdata)
res <- RFLPdist(RFLPdata)
names(res) ## number of bands
str(res$"6")


###################################################
### code chunk number 4: other
###################################################
res1 <- RFLPdist(RFLPdata, distfun = function(x) dist(x, method = "manhattan"))
res2 <- RFLPdist(RFLPdata, distfun = function(x) dist(x, method = "maximum"))
str(res[[1]])
str(res1[[1]])
str(res2[[1]])


###################################################
### code chunk number 5: cor
###################################################
library(MKmisc)
res3 <- RFLPdist(RFLPdata, distfun = corDist)
str(res3$"9")


###################################################
### code chunk number 6: hclust
###################################################
par(mfrow = c(2,2))
plot(hclust(res[[1]]), main = "Euclidean distance")
plot(hclust(res1[[1]]), main = "Manhattan distance")
plot(hclust(res2[[1]]), main = "Maximum distance")
plot(hclust(res3[[1]]), main = "Pearson correlation distance")


###################################################
### code chunk number 7: cutree
###################################################
clust4bd <- hclust(res[[2]])
cgroups50 <- cutree(clust4bd, h=50)
cgroups50


###################################################
### code chunk number 8: sim1
###################################################
library(RColorBrewer)
library(MKmisc)
myCol <- colorRampPalette(brewer.pal(8, "RdYlGn"))(128)
ord <- order.dendrogram(as.dendrogram(hclust(res[[1]])))
temp <- as.matrix(res[[1]])
simPlot(temp[ord,ord], col = rev(myCol), minVal = 0, 
        labels = colnames(temp), title = "(Dis-)Similarity Plot")


###################################################
### code chunk number 9: sim2
###################################################
library(lattice)
print(levelplot(temp[ord,ord], col.regions = rev(myCol),
          at = do.breaks(c(0, max(temp)), 128),
          xlab = "", ylab = "",
          ## Rotate labels of x-axis
          scales = list(x = list(rot = 90)),
          main = "(Dis-)Similarity Plot"))


###################################################
### code chunk number 10: sim3
###################################################
## Euclidean distance
data(RFLPdata)
data(RFLPref)
nrBands(RFLPdata)
res0 <- RFLPdist(RFLPdata, nrBands = 9)
res1 <- RFLPdist2(RFLPdata, nrBands = 9, nrMissing = 1)
res2 <- RFLPdist2(RFLPdata, nrBands = 9, nrMissing = 2)
res3 <- RFLPdist2(RFLPdata, nrBands = 9, nrMissing = 3)


###################################################
### code chunk number 11: sim4
###################################################
par(mfrow = c(2,2))
plot(hclust(res0), main = "0 bands missing")
plot(hclust(res1), main = "1 band missing")
plot(hclust(res2), main = "2 bands missing")
plot(hclust(res3), main = "3 bands missing")


###################################################
### code chunk number 12: RFLPlod
###################################################
RFLPdata.lod <- RFLPlod(RFLPdata, LOD = 60)
par(mfrow = c(1, 2))
RFLPplot(RFLPdata, nrBands = 4, ylim = c(40, 670))
RFLPplot(RFLPdata.lod, nrBands = 4, ylim = c(40, 670))
title(sub = "After applying RFLPlod")


###################################################
### code chunk number 13: sim5
###################################################
res0 <- RFLPdist(RFLPdata, nrBands = 4)
res1.lod <- RFLPdist2(RFLPdata, nrBands = 4, nrMissing = 1, LOD = 60)
ord <- order.dendrogram(as.dendrogram(hclust(res1.lod)))
temp <- as.matrix(res1.lod)
simPlot(temp[ord,ord], col = rev(myCol), minVal = 0, 
        labels = colnames(temp), 
        title = "(Dis-)Similarity Plot\n1 band missing below LOD")


###################################################
### code chunk number 14: RFLPrefplot
###################################################
RFLPrefplot(RFLPdata, RFLPref, nrBands = 9, cex.axis = 0.8)


###################################################
### code chunk number 15: read.blast
###################################################
Dir <- system.file("extdata", package = "RFLPtools") # input directory 
filename <- file.path(Dir, "BLASTexample.txt")
BLAST1 <- read.blast(file = filename)
str(BLAST1)


###################################################
### code chunk number 16: blast
###################################################
data(BLASTdata)


###################################################
### code chunk number 17: simMatrix (eval = FALSE)
###################################################
## res <- simMatrix(BLASTdata)


###################################################
### code chunk number 18: blast1
###################################################
res1 <- simMatrix(BLASTdata, sequence.range = TRUE, Min = 100, Max = 450)
res2 <- simMatrix(BLASTdata, sequence.range = TRUE, Min = 500)


###################################################
### code chunk number 19: blast2
###################################################
library(MKmisc)
simPlot(res2, col = myCol, minVal = 0, cex.axis = 0.5,
        labels = colnames(res2), title = "(Dis-)Similarity Plot")


###################################################
### code chunk number 20: blast3
###################################################
library(lattice)
txt <- trellis.par.get("add.text")
txt$cex <- 0.5
trellis.par.set("add.text" = txt)
myCol <- colorRampPalette(brewer.pal(8, "RdYlGn"))(128)


###################################################
### code chunk number 21: blast31
###################################################
print(levelplot(res2, col.regions = myCol,
          at = do.breaks(c(0, max(res2)), 128),
          xlab = "", ylab = "", 
          ## Rotate labels of x axis
          scales = list(x = list(rot = 90)),
          main = "(Dis-)Similarity Plot"))


###################################################
### code chunk number 22: blast4
###################################################
res.d <- sim2dist(res2)


###################################################
### code chunk number 23: blast5
###################################################
## hierarchical clustering
plot(hclust(res.d), cex = 0.7)


