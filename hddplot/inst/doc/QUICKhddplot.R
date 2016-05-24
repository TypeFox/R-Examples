
## ----setup, cache=FALSE, echo=FALSE--------------------------------------
library(knitr)
options(replace.assign=FALSE,width=72)
opts_chunk$set(fig.path='figs/key-', cache.path='cache/hdd-',
               fig.align='center', dev='pdf', fig.width=3.75,
               fig.height=4, out.width="0.47\\textwidth",
               fig.show='hold', par=TRUE,
               tidy=FALSE,  comment=NA)
knit_hooks$set(par=function(before, options, envir){
if (before && options$fig.show!='none') par(mar=c(4,4,1.6,.1),
              cex.lab=.95,cex.axis=.9,mgp=c(2,.7,0),tcl=-.3)
}, crop=hook_pdfcrop)
pdf.options(pointsize=12)
oldopt <- options(digits=4)


## ----gpsOFinterest-------------------------------------------------------
library(hddplot)
data(golubInfo)
with(golubInfo, table(cancer, tissue.mf))


## ----identify-samples----------------------------------------------------
attach(golubInfo)
## Identify allB samples for that are BM:f or BM:m or PB:m
subsetB <- cancer=="allB" & tissue.mf%in%c("BM:f","BM:m","PB:m")
## Form vector that identifies these as BM:f or BM:m or PB:m
tissue.mfB <- tissue.mf[subsetB, drop=TRUE]
## Separate off the relevant columns of the matrix Golub
data(Golub)          # NB: variables (rows) by cases (columns)
GolubB <- Golub[,subsetB]
detach(golubInfo)


## ----Boxplots, echo=FALSE------------------------------------------------
## Display distributions for the first 20 observations
boxplot(data.frame(GolubB[, 1:20]))  # First 20 columns (observations)
## Random selection of 20 rows (features)
boxplot(data.frame(GolubB[sample(1:7129, 20), ]))


## ----Boxplots, eval=FALSE------------------------------------------------
## ## Display distributions for the first 20 observations
## boxplot(data.frame(GolubB[, 1:20]))  # First 20 columns (observations)
## ## Random selection of 20 rows (features)
## boxplot(data.frame(GolubB[sample(1:7129, 20), ]))


## ----Flawed-scores, eval=TRUE, echo=FALSE--------------------------------
## Uses orderFeatures() (hddplot); see below
ord15 <- orderFeatures(GolubB, cl=tissue.mfB)[1:15]
## Panel A
dfB.ord <- data.frame(t(GolubB[ord15, ]))
## Calculations for the left panel
## Transpose to observations by features
dfB15 <- data.frame(t(GolubB[ord15, ]))
library(MASS)
dfB15.lda <-  lda(dfB15, grouping=tissue.mfB)
scores <- predict(dfB15.lda, dimen=2)$x
## Scores for the single PB:f observation
df.PBf <- with(golubInfo,
  data.frame(t(Golub[ord15, tissue.mf=="PB:f" & cancer=="allB",
                     drop=FALSE])))
scores.PBf <- predict(dfB15.lda, newdata=df.PBf, dimen=2)$x
## For comparison: simulated scores
simscores <- simulateScores(nrow=7129, cl=rep(1:3, c(19,10,2)),
                            cl.other=4, nfeatures=15, seed=41)
  # Returns list elements: scores, cl, scores.other & cl.other


## ----misleading-plots, echo=FALSE----------------------------------------
opar <- par(mar=c(4,4,2.6,.1))
## Warning! The plot that now follows may be misleading!
## Use scoreplot(), from the hddplot package
scoreplot(list(scores=scores, cl=tissue.mfB, other=scores.PBf,
               cl.other="PB:f"))
## Panel B: Repeat plot, now with random normal data
scoreplot(simscores)
par(opar)


## ----Flawed-scores, eval=FALSE, echo=TRUE--------------------------------
## ## Uses orderFeatures() (hddplot); see below
## ord15 <- orderFeatures(GolubB, cl=tissue.mfB)[1:15]
## ## Panel A
## dfB.ord <- data.frame(t(GolubB[ord15, ]))
## ## Calculations for the left panel
## ## Transpose to observations by features
## dfB15 <- data.frame(t(GolubB[ord15, ]))
## library(MASS)
## dfB15.lda <-  lda(dfB15, grouping=tissue.mfB)
## scores <- predict(dfB15.lda, dimen=2)$x
## ## Scores for the single PB:f observation
## df.PBf <- with(golubInfo,
##   data.frame(t(Golub[ord15, tissue.mf=="PB:f" & cancer=="allB",
##                      drop=FALSE])))
## scores.PBf <- predict(dfB15.lda, newdata=df.PBf, dimen=2)$x
## ## For comparison: simulated scores
## simscores <- simulateScores(nrow=7129, cl=rep(1:3, c(19,10,2)),
##                             cl.other=4, nfeatures=15, seed=41)
##   # Returns list elements: scores, cl, scores.other & cl.other


## ----misleading-plots, eval=FALSE----------------------------------------
## opar <- par(mar=c(4,4,2.6,.1))
## ## Warning! The plot that now follows may be misleading!
## ## Use scoreplot(), from the hddplot package
## scoreplot(list(scores=scores, cl=tissue.mfB, other=scores.PBf,
##                cl.other="PB:f"))
## ## Panel B: Repeat plot, now with random normal data
## scoreplot(simscores)
## par(opar)


## ----F-stats, results="hide"---------------------------------------------
## In the following, B is too small for the simulation to give a
## good indication of behaviour in the extreme tail.
library(multtest, quietly=TRUE)
GolubB.maxT <- mt.maxT(GolubB, unclass(tissue.mfB)-1, test="f",
                       B=1000)


## ----qq-Fstats, eval=FALSE-----------------------------------------------
## ## Compare calculated F-statistics with permutation distribution
## qqthin(qf(1-GolubB.maxT$rawp, 2, 28), GolubB.maxT$teststat,
##        print.thinning.details = FALSE)
## ## Compare calculated F-statistics with theoretical F-distribution
## qqthin(qf(ppoints(7129), 2, 28), GolubB.maxT$teststat,
##        print.thinning.details = FALSE)
##   # The theoretical F-distribution gives estimates of quantiles
##   # that are too small
## ## NB also the comparison between the permutation distribution
## ## and the theoretical F:
## qqthin(qf(ppoints(7129), 2, 28), qf(1-GolubB.maxT$rawp, 2, 28),
##        print.thinning.details = FALSE)
##   # qqthin() is a version of qqplot() that thins out points where
##   # overlap is substantial, thus giving smaller graphics files.


## ----plot-Fstats, echo=FALSE, fig.width=3, fig.height=3.25, out.width="0.32\\textwidth"----
## Compare calculated F-statistics with permutation distribution
qqthin(qf(1-GolubB.maxT$rawp, 2, 28), GolubB.maxT$teststat,
       print.thinning.details = FALSE)
## Compare calculated F-statistics with theoretical F-distribution
qqthin(qf(ppoints(7129), 2, 28), GolubB.maxT$teststat,
       print.thinning.details = FALSE)
  # The theoretical F-distribution gives estimates of quantiles
  # that are too small
## NB also the comparison between the permutation distribution
## and the theoretical F:
qqthin(qf(ppoints(7129), 2, 28), qf(1-GolubB.maxT$rawp, 2, 28),
       print.thinning.details = FALSE)
  # qqthin() is a version of qqplot() that thins out points where
  # overlap is substantial, thus giving smaller graphics files.


## ----F-stats, eval=FALSE-------------------------------------------------
## ## In the following, B is too small for the simulation to give a
## ## good indication of behaviour in the extreme tail.
## library(multtest, quietly=TRUE)
## GolubB.maxT <- mt.maxT(GolubB, unclass(tissue.mfB)-1, test="f",
##                        B=1000)


## ----qq-Fstats, eval=FALSE-----------------------------------------------
## ## Compare calculated F-statistics with permutation distribution
## qqthin(qf(1-GolubB.maxT$rawp, 2, 28), GolubB.maxT$teststat,
##        print.thinning.details = FALSE)
## ## Compare calculated F-statistics with theoretical F-distribution
## qqthin(qf(ppoints(7129), 2, 28), GolubB.maxT$teststat,
##        print.thinning.details = FALSE)
##   # The theoretical F-distribution gives estimates of quantiles
##   # that are too small
## ## NB also the comparison between the permutation distribution
## ## and the theoretical F:
## qqthin(qf(ppoints(7129), 2, 28), qf(1-GolubB.maxT$rawp, 2, 28),
##        print.thinning.details = FALSE)
##   # qqthin() is a version of qqplot() that thins out points where
##   # overlap is substantial, thus giving smaller graphics files.


## ----selection-4lda------------------------------------------------------
##              Selection of features that discriminate
## ss 12.3.3: Accuracies and Scores for test data
Golub.BM <- with(golubInfo, Golub[, BM.PB=="BM"])
cancer.BM <- with(golubInfo, cancer[BM.PB=="BM"])
## Now split each of the cancer.BM categories between two subsets
##  Uses divideUp(), from hddplot
gp.id <- divideUp(cancer.BM, nset=2, seed=29)
  # Set seed to allow exact reproduction of the results below
table(gp.id, cancer.BM)


## ----train-test----------------------------------------------------------
accboth <- accTrainTest(x = Golub.BM, cl = cancer.BM,
                        traintest=gp.id, , print.progress=FALSE)


## ----plot-train-test, fig.width=7.5, fig.height=4, out.width="0.97\\textwidth", echo=FALSE----
opar <- par(mar=c(4,4,3.1,.1))
## Use function plotTrainTest() from hddplot
plotTrainTest(x=Golub.BM, nfeatures=c(14,10), cl=cancer.BM, traintest=gp.id)
par(opar)


## ----plot-train-test, eval=FALSE-----------------------------------------
## opar <- par(mar=c(4,4,3.1,.1))
## ## Use function plotTrainTest() from hddplot
## plotTrainTest(x=Golub.BM, nfeatures=c(14,10), cl=cancer.BM, traintest=gp.id)
## par(opar)


## ----match---------------------------------------------------------------
rbind(accboth$sub1.2[1:20],accboth$sub2.1[1:20])
match(accboth$sub1.2[1:20],accboth$sub2.1[1:20])


## ----opt-tissue-mfB-cv, warning=FALSE------------------------------------
##  Cross-validation to determine the optimum number of features
## Accuracy measure will be: tissue.mfB.cv$acc.cv
tissue.mfB.cv <- cvdisc(GolubB, cl=tissue.mfB, nfeatures=1:23,
                        nfold=c(5,1), print.progress=FALSE)
  # 5-fold CV (x1)
## Defective measures will be in acc.resub (resubstitution)
## and acc.sel1 (select features prior to cross-validation)
tissue.mfB.badcv <- defectiveCVdisc(GolubB, cl=tissue.mfB,
                                   foldids=tissue.mfB.cv$folds,
                                   nfeatures=1:23, nfold=c(5,1),
                                   print.progress=FALSE)
## NB: Warning messages have been omitted


## ----tissue-random, warning=FALSE----------------------------------------
## Calculations for random normal data:
set.seed(43)
rGolubB <- matrix(rnorm(prod(dim(GolubB))), nrow=dim(GolubB)[1])
rtissue.mfB.cv <- cvdisc(rGolubB, cl=tissue.mfB, nfeatures=1:23,
                         nfold=c(5,1), print.progress=FALSE)
rtissue.mfB.badcv <- defectiveCVdisc(rGolubB, cl=tissue.mfB,
                                   nfeatures=1:23, nfold=c(5,1),
                                   foldids=rtissue.mfB.cv$folds,
                                   print.progress=FALSE)


## ----plot-acc------------------------------------------------------------
## This function will be used for the plots
plot.acc <- function(cv=cv1, badcv=badcv1, nseq=NULL, badnseq=NULL,
                     title="", ylab="Predictive accuracy",
                     add.legend=TRUE){
  maxg <- min(c(length(badcv$acc.resub), length(cv$acc.cv)))
  if(is.null(nseq))nseq <- 1:maxg
  plot(nseq, badcv$acc.resub[1:maxg], ylim=c(0,1), type="n",
       yaxs="i", xlab="Number of features selected", ylab=ylab)
  par(xpd=T)
  points(nseq, badcv$acc.resub[1:maxg], col=2, type="b", lty=2,
         pch=0, cex=0.8)
  par(xpd=FALSE)
  points(nseq, badcv$acc.sel1[1:maxg], col="gray40", pch=3, cex=0.8)
  lines(lowess(nseq, badcv$acc.sel1[1:maxg], f=.325, iter=0),
        col="gray40", lty=2)
  points(nseq, cv$acc.cv[1:maxg], col="blue", pch=1, cex=0.8)
  lines(lowess(nseq, cv$acc.cv[1:maxg], f=.325, iter=0), col="blue",
        lwd=2)
  xy <- par()$usr[c(1,3)]
  if(add.legend)
    legend(xy[1], xy[2], xjust=0, yjust=0,
           legend=c("Training set 'accuracy'",
             "Defective cross-validation",
             "Cross-validation - select at each fold"),
           lty=c(1,2,1), lwd=c(1,1,2), pch=c(0,3,1),
           col=c("red","gray40","blue"), cex=0.875)
  mtext(side=3,line=0.35, title, adj=0)
}


## ----cv-bad--------------------------------------------------------------
plot.acc(tissue.mfB.cv, tissue.mfB.badcv,
         title="A: Golub data (as for Figure 12.9)")
plot.acc(rtissue.mfB.cv, rtissue.mfB.badcv, ylab="",
         title="B: Random data", add.legend=FALSE)


## ----which---------------------------------------------------------------
##                          Which features?
genelist <- matrix(tissue.mfB.cv$genelist[1:3, ,], nrow=3)
tab <- table(genelist, row(genelist))
ord <- order(tab[,1], tab[,2], decreasing=TRUE)
tab[ord,]


## ----cv-BMsamples--------------------------------------------------------
##         Cross-validation: bone marrow ({BM}) samples only
BMonly.cv <- cvdisc(Golub.BM, cl=cancer.BM, nfeatures=1:25,
                    nfold=c(5,1), print.progress=FALSE)
tissue.mfB.scores <-
  cvscores(cvlist = tissue.mfB.cv, nfeatures = 3, cl.other = NULL,
           print.progress=FALSE)
BMonly.scores <- cvscores(cvlist=BMonly.cv, nfeatures=19,
                          cl.other=NULL, print.progress=FALSE)


## ----cv-Bcell-gphAB, echo=FALSE------------------------------------------
opar <- par(mar=c(4,4,2.6,.1))
## Panel A: Uses tissue.mfB.acc from above
scoreplot(scorelist = tissue.mfB.scores, cl.circle=NULL,
          prefix="B-cell subset -")
## Panel B; classify bone marrow samples a/c cancer type.
scoreplot(scorelist=BMonly.scores, cl.circle=tissue.mfB,
          circle=tissue.mfB%in%c("BM:f","BM:m"),
          params=list(circle=list(col=c("cyan","gray"))),
          prefix="B: BM samples -")
par(opar)


## ----cv-Bcell-gphAB, eval=FALSE------------------------------------------
## opar <- par(mar=c(4,4,2.6,.1))
## ## Panel A: Uses tissue.mfB.acc from above
## scoreplot(scorelist = tissue.mfB.scores, cl.circle=NULL,
##           prefix="B-cell subset -")
## ## Panel B; classify bone marrow samples a/c cancer type.
## scoreplot(scorelist=BMonly.scores, cl.circle=tissue.mfB,
##           circle=tissue.mfB%in%c("BM:f","BM:m"),
##           params=list(circle=list(col=c("cyan","gray"))),
##           prefix="B: BM samples -")
## par(opar)


