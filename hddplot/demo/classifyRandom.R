require(MASS)
require(hddplot)

## Generate random data
gps <- rep(1:3, c(19,10,5))
ExpMat <- matrix(rnorm(7129*length(gps)), nrow=7129)
simscores <- simulateScores(x=ExpMat, cl=gps, nfeatures=15, dimen=2)
scoreplot(simscores, prefix.title="Random normal data:")

## The 15 features (from 7129) were chosen that individually best discriminated.
## The discriminant scores are for the data used to derive the scores.
##

## Now vary number of features that are used

goodCV <- cvdisc(ExpMat, cl=gps, nfeatures=1:23, nfold=c(5,2))
badCV <- defectiveCVdisc(ExpMat, cl=gps, nfeatures=1:23, nfold=c(5,2))
maxg <- min(c(length(badCV$acc.resub), length(goodCV$acc.cv)))
nseq <- 1:maxg
plot(nseq, badCV$acc.resub[1:maxg], ylim=c(0,1), type="n", yaxs="i",
     xlab="Number of features selected", ylab="Predictive 'accuracy'")
par(xpd=T)
points(nseq, badCV$acc.resub[1:maxg], col=2, type="b", lty=2, pch=0,
       cex=0.8)
par(xpd=FALSE)
points(nseq, badCV$acc.sel1[1:maxg], col="gray40", pch=3, cex=0.8)
lines(lowess(nseq, badCV$acc.sel1[1:maxg], f=.325, iter=0), col="gray40",
      lty=2)
points(nseq, goodCV$acc.cv[1:maxg], col="blue", pch=1, cex=0.8)
lines(lowess(nseq, goodCV$acc.cv[1:maxg], f=.325, iter=0), col="blue",lwd=2)
xy <- par()$usr[c(1,3)]
legend(xy[1], xy[2], xjust=0, yjust=0,
       legend=c("Training set accuracy", "Defective cross-validation",
                "Cross-validation - select at each fold"),
       lty=c(1,2,1), lwd=c(1,1,2), pch=c(0,3,1),
       col=c("red","gray40","blue"), cex=0.875)
title <- "Random data; 'accuracy' measures vs # of features selected"
mtext(side=3,line=0.35, title, adj=0)

## Training set accuracy uses the training data to 'estimate' accuracy.
## Defective CV uses all the data to select the 'best' features.
## Cross-validation repeats the selection process at each fold.
