
## ----preliminary, echo=FALSE, results="hide", message=FALSE, fig.keep="none"----
options(width=60, prompt="R> ", continue="     ", useFancyQuotes=FALSE, digits=3)
library(WeightedCluster)
library(TraMineR)
library(knitr)
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
## knit_hooks$set(crop =hook_pdfcrop)
knit_hooks$set(small.mar = function(before, options, envir) {
    if (before)  par(mar=c(2.1, 4.1, 4.1, 1.1))  # smaller margin on top and right
})
opts_chunk$set(message=FALSE, prompt=TRUE, dev="pdf", echo=TRUE, comment=NA, small.mar=TRUE, fig.align="center", fig.path="graphics/WCP-", size="small")
## knit_hooks$set(error = function(x, options) stop(x))
##knit_hooks$set(warning = function(x, options) stop("Warnings: ", x))


## ----install, echo=TRUE, results="hide", eval=FALSE, fig.keep="none"-----
## install.packages("WeightedCluster")
## library(WeightedCluster)


## ----dataload, echo=TRUE, results="hide", fig.keep="none"----------------
data(mvad)


## ----wcagg, echo=TRUE, fig.keep="none"-----------------------------------
aggMvad <- wcAggregateCases(mvad[, 17:86])
print(aggMvad)
uniqueMvad <- mvad[aggMvad$aggIndex, 17:86]


## ----wcagg-diss, echo=TRUE, fig.keep="none"------------------------------
mvad.seq <- seqdef(uniqueMvad, weights=aggMvad$aggWeights)
## Computing Hamming distance between sequence
diss <- seqdist(mvad.seq, method="HAM")


## ----hierclust, echo=TRUE, consolew=50, fig.keep="none"------------------
averageClust <- hclust(as.dist(diss), method="average", members=aggMvad$aggWeights)


## ----avgtreecomputeecho, echo=TRUE, eval=FALSE, fig.keep="none"----------
## averageTree <- as.seqtree(averageClust, seqdata=mvad.seq, diss=diss, ncluster=6)
## seqtreedisplay(averageTree, type="d", border=NA,  showdepth=TRUE)


## ----avgqualcompute, echo=TRUE, consolew=40, fig.keep="none"-------------
avgClustQual <- as.clustrange(averageClust, diss, weights=aggMvad$aggWeights, ncluster=10)


## ----avgqualplot, echo=TRUE, fig.height=3, fig.width=9-------------------
plot(avgClustQual)


## ----avgqualplotnorm, echo=TRUE, fig.height=3, fig.width=9---------------
plot(avgClustQual, norm="zscore")


## ----avgqualprint, echo=TRUE, fig.keep="none"----------------------------
summary(avgClustQual, max.rank=2)


## ----pamcompute, echo=TRUE, consolew=40, fig.keep="none"-----------------
pamClustRange <- wcKMedRange(diss, kvals=2:10, weights=aggMvad$aggWeights)


## ----pamqualprint, echo=TRUE, fig.keep="none"----------------------------
summary(pamClustRange, max.rank=2)


## ----mppamplot, fig.height=6, fig.width=9--------------------------------
seqdplot(mvad.seq, group=pamClustRange$clustering$cluster5, border=NA)


## ----mvadcluster, echo=TRUE, fig.keep="none"-----------------------------
uniqueCluster5 <- avgClustQual$clustering$cluster5
mvad$cluster5 <- uniqueCluster5[aggMvad$disaggIndex]
table(mvad$funemp, mvad$cluster5)


