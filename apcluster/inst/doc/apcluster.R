## ----Init,echo=FALSE,message=FALSE,results='hide'---------------------
options(width=72)
knitr::opts_knit$set(width=72)
set.seed(0)
library(kebabs, quietly=TRUE)
library(apcluster, quietly=TRUE)
apclusterVersion <- packageDescription("apcluster")$Version
apclusterDateRaw <- packageDescription("apcluster")$Date
apclusterDateYear <- as.numeric(substr(apclusterDateRaw, 1, 4))
apclusterDateMonth <- as.numeric(substr(apclusterDateRaw, 6, 7))
apclusterDateDay <- as.numeric(substr(apclusterDateRaw, 9, 10))
apclusterDate <- paste(month.name[apclusterDateMonth], " ",
                     apclusterDateDay, ", ",
                     apclusterDateYear, sep="")

## ----InstallAPCluster,eval=FALSE--------------------------------------
#  install.packages("apcluster")

## ----LoadAPCluster,eval=FALSE-----------------------------------------
#  library(apcluster)

## ----OpenVignette,eval=FALSE------------------------------------------
#  vignette("apcluster")

## ----ShowHelp,eval=FALSE----------------------------------------------
#  help(apcluster)

## ----CreateDataSet1,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
cl1 <- cbind(rnorm(30, 0.3, 0.05), rnorm(30, 0.7, 0.04))
cl2 <- cbind(rnorm(30, 0.7, 0.04), rnorm(30, 0.4, .05))
x1 <- rbind(cl1, cl2)
plot(x1, xlab="", ylab="", pch=19, cex=0.8)

## ----APClusterDataSet1------------------------------------------------
apres1a <- apcluster(negDistMat(r=2), x1)

## ----APClusterDataSet1b-----------------------------------------------
s1 <- negDistMat(x1, r=2)
apres1b <- apcluster(s1)

## ----ShowHelpAPResult,eval=FALSE--------------------------------------
#  help(APResult)

## ----ShowResultAPClusterDataSet1--------------------------------------
apres1a

## ----PlotResultAPClusterDataSet1,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
plot(apres1a, x1)

## ----HeatmapResultAPClusterDataSet1,fig.width=7,fig.height=7,out.width='0.6\\textwidth'----
heatmap(apres1a)

## ----HeatmapResultAPClusterDataSet1b,fig.width=7,fig.height=7,out.width='0.6\\textwidth'----
heatmap(apres1b, s1)

## ----APClusterDataSet1Details-----------------------------------------
apres1c <- apcluster(s1, details=TRUE)

## ----PlotAPClusterDataSet1Details,fig.width=6,fig.height=4,out.width='0.5\\textwidth'----
plot(apres1c)

## ----APClusterDataSet1convits15---------------------------------------
apres1c <- apcluster(s1, convits=15, details=TRUE)
apres1c

## ----CreateDataSet2,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
cl3 <- cbind(rnorm(20, 0.50, 0.03), rnorm(20, 0.72, 0.03))
cl4 <- cbind(rnorm(25, 0.50, 0.03), rnorm(25, 0.42, 0.04))
x2 <- rbind(x1, cl3, cl4)
plot(x2, xlab="", ylab="", pch=19, cex=0.8)

## ----APClusterDataSet2,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
apres2a <- apcluster(negDistMat(r=2), x2)
plot(apres2a, x2)

## ----APClusterDataSet2q0,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
apres2b <- apcluster(negDistMat(r=2), x2, q=0)
plot(apres2b, x2)

## ----PlotAPClusterDataSet2q08,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
apres2c <- apcluster(negDistMat(r=2), x2, q=0.8)
plot(apres2c, x2)

## ----APClusterDataSet2q08showp----------------------------------------
apres2c@p

## ----HeatmapResultAPClusterDataSet2q08,fig.width=7,fig.height=7,out.width='0.6\\textwidth'----
heatmap(apres2c)

## ----PreferenceRangeDataSet2------------------------------------------
preferenceRange(apres2b@sim)

## ----APClusterKDataSet2,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
apres2d <- apclusterK(negDistMat(r=2), x2, K=2, verbose=TRUE)
plot(apres2d, x2)

## ----IrisData1--------------------------------------------------------
data(iris)
apIris1 <- apcluster(negDistMat(r=2), iris)
apIris1

## ----IrisDataPlot1,fig.width=10,fig.height=10,out.width='\\textwidth'----
plot(apIris1, iris)

## ----IrisDataHeatmap1,fig.width=7,fig.height=7,out.width='0.6\\textwidth'----
heatmap(apIris1)

## ----IrisData2--------------------------------------------------------
data(iris)
apIris2 <- apcluster(negDistMat(r=2), iris, q=0)
apIris2

## ----IrisDataPlot,fig.width=10,fig.height=10,out.width='\\textwidth'----
plot(apIris2, iris)

## ----IrisDataHeatmap2,fig.width=7,fig.height=7,out.width='0.6\\textwidth'----
heatmap(apIris2)

## ----AggExClusterDataSet1---------------------------------------------
aggres1a <- aggExCluster(negDistMat(r=2), x1)
aggres1a

## ----DendrogramAggExClusterDataSet1,fig.width=5,fig.height=5,out.width='0.5\\textwidth'----
plot(aggres1a)

## ----HeatmapAggExClusterDataSet1,fig.width=7,fig.height=7,out.width='0.6\\textwidth'----
heatmap(aggres1a, s1)

## ----ExtractAggExClustersDataSet1,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
cl1a <- cutree(aggres1a, k=2)
cl1a
plot(cl1a, x1)

## ----AggExClusterAPDataSet2q08----------------------------------------
aggres2a <- aggExCluster(x=apres2c)
aggres2a

## ----DendrogramAggExAPDataSet2,fig.width=5,fig.height=5,out.width='0.5\\textwidth'----
plot(aggres2a)

## ----DendrogramAggExAPDataSet2b,fig.width=5,fig.height=5,out.width='0.5\\textwidth'----
plot(aggres2a, showSamples=TRUE, nodePar=list(pch=NA, lab.cex=0.4))

## ----HeatmapAggExAPDataSet2,fig.width=7,fig.height=7,out.width='0.6\\textwidth'----
heatmap(aggres2a)

## ----PlotAggExAPDataSet2k25,fig.width=8,fig.height=8,out.width='\\textwidth'----
par(mfrow=c(2,2))
for (k in 5:2)
    plot(aggres2a, x2, k=k, main=paste(k, "clusters"))

## ----APClusterLevDataSet3,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
cl5 <- cbind(rnorm(100, 0.3, 0.05), rnorm(100, 0.7, 0.04))
cl6 <- cbind(rnorm(100, 0.70, 0.04), rnorm(100, 0.4, 0.05))
x3 <- rbind(cl5, cl6)
apres3 <- apclusterL(s=negDistMat(r=2), x=x3, frac=0.1, sweeps=5, p=-0.2)
apres3
plot(apres3, x3)

## ----APClusterLevResultDataSet3---------------------------------------
dim(apres3@sim)
apres3@sel
apres3@netsimLev

## ----APClusterLevDataSet3Heat,fig.width=7,fig.height=7,out.width='0.6\\textwidth'----
heatmap(apres3)

## ----SparseEx1--------------------------------------------------------
dsim <- negDistMat(x2, r=2)
ssim <- as.SparseSimilarityMatrix(dsim, lower=-0.2)
str(ssim)

## ----SparseEx1Run,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
sapres <- apcluster(ssim, q=0)
plot(sapres, x2)

## ----SparseEx1Run2----------------------------------------------------
preferenceRange(ssim)
apclusterK(ssim, K=2)

## ----SparseEx1RunHeatmap,fig.width=7,fig.height=7,out.width='0.6\\textwidth'----
heatmap(sapres, ssim)

## ----LoadCh22Promoters------------------------------------------------
library(Biostrings)
filepath <- system.file("examples", "ch22Promoters.fasta",
                        package="apcluster")
ch22Promoters <- readDNAStringSet(filepath)
ch22Promoters

## ----SimCh22Promoters-------------------------------------------------
library(kebabs)
specK6 <- spectrumKernel(k=6)
promSim <- getKernelMatrix(specK6, ch22Promoters)
promSim[1:4, 1:4]

## ----APCh22Promoters--------------------------------------------------
promAP <- apcluster(promSim, q=0)
promAP

## ----HeatMapAPCh22Promoters,fig.width=7,fig.height=7,out.width='0.6\\textwidth'----
heatmap(promAP, promSim, Rowv=FALSE, Colv=FALSE)

## ----aggExCh22Promoters-----------------------------------------------
promAgg <- aggExCluster(promSim, promAP)

## ----DendrogramAPCh22Promoters,fig.width=5,fig.height=5,out.width='0.5\\textwidth'----
plot(promAgg)

## ----ExtractAggCh22Promoters------------------------------------------
prom5 <- cutree(promAgg, k=5)
prom5

## ----HeatMap5Ch22Promoters,fig.width=7,fig.height=7,out.width='0.6\\textwidth'----
heatmap(prom5, promSim, Rowv=FALSE, Colv=FALSE)

## ----NegDistMatDataSet2-----------------------------------------------
s <- negDistMat(x2)

## ----CreateToyData----------------------------------------------------
ex <- matrix(c(0, 0.5, 0.8, 1, 0, 0.2, 0.5, 0.7,
               0.1, 0, 1, 0.3, 1, 0.8, 0.2), 5, 3,byrow=TRUE)
ex

## ----NegEuclDistMatToyData--------------------------------------------
negDistMat(ex)

## ----NegSqEuclDistMatToyData------------------------------------------
negDistMat(ex, r=2)

## ----NegMaxDistToyData------------------------------------------------
negDistMat(ex, method="maximum")

## ----NegManhattanDistToyData------------------------------------------
negDistMat(ex, method="manhattan")

## ----NegCanberraDistToyData-------------------------------------------
negDistMat(ex, method="canberra")

## ----NegMinkowskiDistToyData------------------------------------------
negDistMat(ex, method="minkowski", p=3)

## ----GetFunction------------------------------------------------------
sim <- negDistMat(r=2)
is.function(sim)
apcluster(sim, x1)

## ----RBFKernelToyData-------------------------------------------------
expSimMat(ex)

## ----LaplaceKernelToyData---------------------------------------------
expSimMat(ex, r=1)

## ----PearsonToyData---------------------------------------------------
corSimMat(ex, method="pearson")

## ----SpearmanToyData--------------------------------------------------
corSimMat(ex, method="spearman")

## ----TruncDistToyData-------------------------------------------------
linSimMat(ex, w=1.2)

## ----LinKernelToyData-------------------------------------------------
linKernel(ex[2:5,])

## ----NormLinKernelToyData---------------------------------------------
linKernel(ex[2:5,], normalize=TRUE)

## ----RectangularNegDistMatDataSet1------------------------------------
sel <- sort(sample(1:nrow(x1), ceiling(0.08 * nrow(x1))))
sel
s1r <- negDistMat(x1, sel, r=2)
dim(s1r)
s1r[1:7,]

## ----SimFunCh22Promoters----------------------------------------------
spectrumK6 <- function(x, sel=NA)
{
    if (any(is.na(sel)))
        s <- getKernelMatrix(specK6, x)
    else
        s <- getKernelMatrix(specK6, x, x[sel])

    as(s, "matrix")
}

## ----APCh22Promoters2-------------------------------------------------
promAPL <- apclusterL(s=spectrumK6, ch22Promoters, frac=0.1, sweeps=10,
                      p=promAP@p)
promAPL

## ----CustomSimSparse--------------------------------------------------
sparseSim <- function(x)
{
    as.SparseSimilarityMatrix(negDistMat(x, r=2), lower=-0.2)
}

sapres2 <- apcluster(sparseSim, x2, q=0)
sapres2
str(similarity(sapres2))

## ----CreateLabeledToyData---------------------------------------------
x3 <- c(1, 2, 3, 7, 8, 9)
names(x3) <- c("a", "b", "c", "d", "e", "f")
s3 <- negDistMat(x3, r=2)

## ----ShowToyDataLabels------------------------------------------------
s3
colnames(s3)

## ----ClusterLabeledToyData--------------------------------------------
apres3a <-apcluster(s3)
apres3a
apres3a@exemplars
apres3a@clusters

## ----ExtractLabelsFromClusterToyData----------------------------------
apres3a@exemplars
labels(apres3a, type="names")
labels(apres3a, type="exemplars")
labels(apres3a, type="enum")

## ----HeatmapResultAPClusterDataSetq08b,fig.width=7,fig.height=7,out.width='0.6\\textwidth'----
heatmap(apres2c, sideColors=c("darkgreen", "yellowgreen"),
        col=terrain.colors(12), Rowv=FALSE, dendScale=0.5,
        margins=c(3, 3, 2), legend="col")

## ----HeatmapResultAPClusterDataSet2q08c,fig.width=7,fig.height=7,out.width='0.6\\textwidth'----
heatmap(apres2c, sideColors=rainbow(length(apres2c)), Rowv=FALSE, Colv=FALSE,
        cexRow=(0.2 + 1 / log10(nrow(apres2c@sim))),
        cexCol=(0.2 + 1 / log10(nrow(apres2c@sim))))

## ----PlotAddLegend,fig.width=5,fig.height=5.5,out.width='0.5\\textwidth'----
plot(apres2a, x2)
legend("bottomleft", legend=paste("Cluster", 1:length(apres2a)),
       col=rainbow(length(apres2a)), pch=19)

## ----PlotOnlyLegend,fig.width=5,fig.height=2.5,out.width='0.5\\textwidth'----
plot.new()
par(oma=rep(0, 4), mar=rep(0, 4))
legend("center", legend=paste("Cluster", 1:length(apres2c)),
       col=rainbow(length(apres2c)), pch=19)

## ----GetBibTeX,eval=FALSE---------------------------------------------
#  toBibtex(citation("apcluster"))

