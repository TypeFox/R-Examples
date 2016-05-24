
## ----preliminary, echo=FALSE, results="hide", message=FALSE, warning=FALSE, fig.keep="none"----
options(width=60, prompt="R> ", continue="     ", useFancyQuotes=FALSE, digits=3)
library(WeightedCluster)
library(TraMineR)
library(vegan)
library(knitr)
library(isotone)
hook_setwidth <- local({
    default.width <- 0
	function(before, options, envir){
		if(before) {
            default.width <<- getOption("width")
			options(width = options$consolew)
		} else{
			options(width = default.width)
		}
		return(NULL)
	}
})
knit_hooks$set(consolew =hook_setwidth)
##knit_hooks$set(crop = hook_pdfcrop)
knit_hooks$set(small.mar = function(before, options, envir) {
    if (before)  par(mar=c(2.1, 4.1, 4.1, 1.1))  # smaller margin on top and right
})
opts_knit$set(concordance=TRUE)
opts_chunk$set(message=FALSE, prompt=TRUE, echo=TRUE, dev="pdf", comment=NA, small.mar=TRUE, fig.align="center", fig.path="graphics/WC-", out.width=".8\\linewidth", size="small")
## knit_hooks$set(error = function(x, options) stop(x))
## knit_hooks$set(warning = function(x, options) stop("Warnings: ", x))


## ----WCload, eval=FALSE, fig.keep="none"---------------------------------
## install.packages("WeightedCluster")
## library(WeightedCluster)


## ----distcompute, tidy=FALSE, fig.keep="none"----------------------------
data(mvad)
mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school",
    "training")
mvad.labels <- c("Employment", "Further Education", "Higher Education",
    "Joblessness", "School", "Training")
mvad.scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
mvadseq <- seqdef(mvad[, 17:86], alphabet = mvad.alphabet, states = mvad.scodes,
    labels = mvad.labels, weights=mvad$weight, xtstep=6)
## Defining the custom cost matrix
subm.custom <- matrix(
      c(0, 1, 1, 2, 1, 1,
        1, 0, 1, 2, 1, 2,
        1, 1, 0, 3, 1, 2,
        2, 2, 3, 0, 3, 1,
        1, 1, 1, 3, 0, 2,
        1, 2, 2, 1, 2, 0),
      nrow = 6, ncol = 6, byrow = TRUE)
## Computing the OM dissimilarities
mvaddist <- seqdist(mvadseq, method="OM", indel=1.5, sm=subm.custom)


## ----hclustcompute, fig.keep="none"--------------------------------------
wardCluster <- hclust(as.dist(mvaddist), method="ward", members=mvad$weight)


## ----as.seqtreecompute, fig.keep="none"----------------------------------
wardTree <- as.seqtree(wardCluster, seqdata=mvadseq, diss=mvaddist, ncluster=6)


## ----seqtreedisplay, echo=FALSE, results="hide", eval=FALSE, fig.keep="none"----
## seqtreedisplay(wardTree, type="d", border=NA, filename="wardtree.png", showdepth=TRUE, showtree=FALSE)


## ----seqtreedisplay-fake, eval=FALSE, echo=TRUE, results="hide", fig.keep="none"----
## seqtreedisplay(wardTree, type="d", border=NA, showdepth=TRUE)


## ----cutreecompute, fig.keep="none"--------------------------------------
clust4 <- cutree(wardCluster, k=4)


## ----seqdplot-4clust, fig.width=7, fig.height=5--------------------------
seqdplot(mvadseq, group=clust4, border=NA)


## ----wcKMedoids-compute, fig.keep="none"---------------------------------
pamclust4 <- wcKMedoids(mvaddist, k=4, weights=mvad$weight)


## ----pamclust4-plot, fig.width=7, fig.height=5---------------------------
seqdplot(mvadseq, group=pamclust4$clustering, border=NA)


## ----wcKMedoids-print, fig.keep="none"-----------------------------------
print(mvadseq[unique(pamclust4$clustering), ], format="SPS")


## ----wcKMedoids-compute4-ward, fig.keep="none"---------------------------
pamwardclust4 <- wcKMedoids(mvaddist, k=4, weights=mvad$weight, initialclust=wardCluster)


## ----pamwardclust4-plot-ward, fig.width=7, fig.height=5------------------
seqdplot(mvadseq, group=pamwardclust4$clustering, border=NA)


## ----wcClusterQuality-compute, echo=2:3, consolew=70, fig.keep="none"----
options(digits=2)
clustqual4 <- wcClusterQuality(mvaddist, clust4, weights=mvad$weight)
clustqual4$stats
options(digits=3)


## ----wcClusterQuality-computeasw, fig.keep="none"------------------------
clustqual4$ASW


## ----silhouette-indexplot, echo=2:3, dev="png", dpi=600, fig.width=7, fig.height=5----
par(mar=c(2.1, 4.1, 4.1, 1.1))
sil <- wcSilhouetteObs(mvaddist, clust4, weights=mvad$weight, measure="ASWw")
seqIplot(mvadseq, group=clust4, sortv=sil)


## ----wcClusterQuality-compute9-clustqual4pam, echo=2, fig.keep="none"----
options(digits=1)
pamclust4$stats
options(digits=3)


## ----wcRange-compute, fig.keep="none"------------------------------------
wardRange <- as.clustrange(wardCluster, diss=mvaddist, weights=mvad$weight, ncluster=20)
summary(wardRange, max.rank=2)


## ----wcRange-plot, fig.width=8, fig.height=3.5---------------------------
plot(wardRange, stat=c("ASWw", "HG", "PBC", "HC"))


## ----wcRange-plotzscore, fig.width=8, fig.height=3.5---------------------
plot(wardRange, stat=c("ASWw", "HG", "PBC", "HC"), norm="zscore")


## ----wardClust6-plot, fig.width=7, fig.height=6--------------------------
seqdplot(mvadseq, group=wardRange$clustering$cluster6, border=NA)


## ----wcKMedRange-compute, echo=TRUE, fig.keep="none"---------------------
pamRange <- wcKMedRange(mvaddist, kvals=2:20, weights=mvad$weight)
summary(pamRange, max.rank=2)


## ----cluster-naming, echo=TRUE, fig.keep="none"--------------------------
mvad$pam4 <- factor(pamclust4$clustering, levels=c(66, 467, 607, 641), labels=c("Appr.-Empl.", "Ecole-Empl.", "Ecole sup.", "Sans empl."))


## ----auto-cluster-naming, echo=TRUE, consolew=80, fig.keep="none"--------
mvad$pam4.auto <- seqclustname(mvadseq, pamclust4$clustering, mvaddist)
table( mvad$pam4.auto, mvad$pam4)


## ----mds-compute, echo=FALSE, fig.keep="none"----------------------------
worsq <- wcmdscale(mvaddist, w=mvad$weight, k=2)
mvad$test <- rep(-1, nrow(mvad))
for(clust in unique(pamclust4$clustering)){
    cond <- pamclust4$clustering == clust
    values <- worsq[cond, 2]
    mvad$test[cond] <- as.integer(values > weighted.median(values, w=mvad$weight[cond]))
}
mvad$test <- factor(mvad$test, levels=0:1, labels=c("non-test", "test"))


## ----testdplot, fig.width=10, fig.height=4.5-----------------------------
seqdplot(mvadseq, group=mvad$test, border=NA)


## ----testiplot-chisq, fig.keep="none"------------------------------------
tb <- xtabs(weight~test+pam4, data=mvad)
chisq.test(tb)


## ----testdplot-EE, fig.width=10, fig.height=4.5--------------------------
EcoleEmploi <- mvad$pam4=="Ecole-Empl."
seqdplot(mvadseq[EcoleEmploi, ], group=mvad$test[EcoleEmploi], border=NA)


## ----testiplot-dissassoc, fig.keep="none"--------------------------------
set.seed(1)
dsa <- dissassoc(mvaddist, mvad$test, weights=mvad$weight, weight.permutation="diss", R=5000)
print(dsa$stat)


## ----seqtree-link-covar, echo=FALSE, results="hide", eval=FALSE, fig.keep="none"----
## set.seed(1)
## tree <- seqtree(mvadseq~gcse5eq+Grammar+funemp, data=mvad, diss=mvaddist, weight.permutation="diss")
## seqtreedisplay(tree, type="d", border=NA, filename="seqtree.png", showtree=FALSE)


## ----seqtree-link-covar-fake, eval=FALSE, echo=TRUE, results="hide", fig.keep="none"----
## tree <- seqtree(mvadseq~gcse5eq+Grammar+funemp, data=mvad, diss=mvaddist, weight.permutation="diss")
## seqtreedisplay(tree, type="d", border=NA)


## ----wcAggregateCases, echo=TRUE, fig.keep="none"------------------------
ac <- wcAggregateCases(mvad[, 17:86], weights=mvad$weight)
ac


## ----wcAggregateCases-seqdef, echo=TRUE, fig.keep="none"-----------------
uniqueSeq <- seqdef(mvad[ac$aggIndex, 17:86], alphabet = mvad.alphabet,
    states = mvad.scodes, labels = mvad.labels,  weights=ac$aggWeights)


## ----wcAggregateCases-wckmedoids, echo=TRUE, fig.keep="none"-------------
mvaddist2 <- seqdist(uniqueSeq, method="OM", indel=1.5, sm=subm.custom)
pamclust4ac <- wcKMedoids(mvaddist2, k=4, weights=ac$aggWeights)


## ----wcAggregateCases-wckmedoids-back, echo=TRUE, fig.keep="none"--------
mvad$acpam4 <- pamclust4ac$clustering[ac$disaggIndex]


## ----wcAggregateCases-wckmedoids-table, echo=TRUE, fig.keep="none"-------
table(mvad$pam4, mvad$acpam4)


## ----perf-loadsimul, echo=FALSE, results="hide", message=FALSE, fig.keep="none"----
load(file="randB.RData")
randB$nnn <- factor(randB$n, levels=c(200, 500, 1000, 2000), labels=c("n=200", "n=500", "n=1000", "n=2000"))
library(lattice)


## ----perf-nbyk-plot-time, fig.width=8, fig.height=3.5, echo=FALSE--------
print(xyplot(ClusterTime~k|nnn, data=randB, groups=test, type="l", auto.key = list(points=F, lines=T,corner = c(0, 0.8), cex=.9, size=2), scales=list(y = list(log = TRUE)), ylab="Temps de calcul", xlab="Nombre de groupes"))


## ----perf-nbyk-plot, fig.width=8, fig.height=3.5, echo=FALSE-------------
print(xyplot(relative.tot~k|nnn, data=randB, groups=test, type="l", auto.key = list(points=F, lines=T,corner = c(0, 0.8), size=2, cex=.9), ylab="Temps relatif", xlab="Nombre de groupes"))


## ----mdscomputeannexe, echo=TRUE, results="hide", fig.keep="none"--------
library(vegan)
worsq <- wcmdscale(mvaddist, w=mvad$weight, k=2)


## ----mdscomputeannexesecond, echo=TRUE, results="hide", fig.keep="none"----
library(isotone)
mvad$test <- rep(-1, nrow(mvad))
for(clust in unique(pamclust4$clustering)){
    cond <- pamclust4$clustering == clust
    values <- worsq[cond, 2]
    mvad$test[cond] <- values> weighted.median(values, w=mvad$weight[cond])
}
mvad$test <- factor(mvad$test, levels=0:1, labels=c("non-test", "test"))


