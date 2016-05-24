### R code from vignette source 'popKorn.Rnw'

###################################################
### code chunk number 1: generate_X
###################################################
p <- 10; n <- 10
set.seed(18)
Xmat <- matrix(rnorm(p*n), nrow=n, ncol=p)
colnames(Xmat) <- paste("p.", 1:p, sep="")


###################################################
### code chunk number 2: plot_X
###################################################
Xmat.melt <- reshape(as.data.frame(Xmat), direction="long", varying=1:p)
colnames(Xmat.melt) <- c("p", "X", "n")
boxplot(X ~p, data=Xmat.melt, col='sienna2', main='Boxplot of sample X',
xlab='Population p', ylab='X')
points(1:10, colMeans(Xmat), pch=20, col='white', cex=2.5)
Xmat.melt$p <- as.factor(Xmat.melt$p)
#sort(colMeans(Xmat), decreasing=TRUE)


###################################################
### code chunk number 3: bonf_1
###################################################
library(popKorn)
bonferroniIntervals(Xmat, k=1)


###################################################
### code chunk number 4: asym_1
###################################################
asymmetricIntervals(Xmat, k=1, eps=0.05)


###################################################
### code chunk number 5: plot_intervals
###################################################
asym.out <- asymmetricIntervals(Xmat, k=3)
bonf.out <- bonferroniIntervals(Xmat, k=3)
library(plotrix)
plotCI(1:3, bonf.out[,1], ui=bonf.out[,"upr"], li=bonf.out[,"lwr"], 
  col="red", scol="darkgreen", slty=2, lwd=2, pch=20, cex=1.8, xaxt='n', 
  xlab="", ylab="")
plotCI(1:3, add=TRUE, asym.out[,1], ui=asym.out[,"upr"], li=asym.out[,"lwr"], 
  scol="darkorange1", lwd=2, pch=NA, sfrac=0.02, gap=0.02)
title("Bonferroni vs. Asymmetric Intervals")
axis(side=1, at=1:3, labels=rownames(asym.out))


###################################################
### code chunk number 6: opt_lam
###################################################
optimalLambdaC(alpha=0.05, n=10, p=10, k=3, var.known=FALSE)


###################################################
### code chunk number 7: generate_Y
###################################################
p <- 10; n <- 10
set.seed(18)
Ymat <- matrix(c(rnorm(n, 3), rnorm((p-1)*n)), nrow=n, ncol=p)
colnames(Ymat) <- paste("p.", 1:p, sep="")


###################################################
### code chunk number 8: plot_Y
###################################################
Ymat.melt <- reshape(as.data.frame(Ymat), direction="long", varying=1:p)
colnames(Ymat.melt) <- c("p", "Y", "n")
boxplot(Y ~ p, data=Ymat.melt, col='sienna2', main='Boxplot of sample Y',
xlab='Population p', ylab='Y')
points(1:10, colMeans(Ymat), pch=20, col='white', cex=2.5)
Ymat.melt$p <- as.factor(Ymat.melt$p)


###################################################
### code chunk number 9: asym_2
###################################################
asymmetricIntervals(Ymat, k=1, eps=0.05)


###################################################
### code chunk number 10: get_theta_diff
###################################################
sample.means <- colMeans(Ymat)
sorted.sample.means <- sort(sample.means, dec=TRUE)
theta.diff <- sorted.sample.means[1:(p-1)] - sorted.sample.means[2:p]


###################################################
### code chunk number 11: true_cov_prob
###################################################
opt.lam.c <- optimalLambdaC(0.05, 10, 10, 1, FALSE, eps=0.05)
exactCoverageProb(theta.diff=theta.diff, c.val=opt.lam.c$c.val, 
lambda=opt.lam.c$lambda, sigma.2=1, n=10)


###################################################
### code chunk number 12: search_for_c
###################################################
grid.c <- seq(0.8*opt.lam.c$c.val, opt.lam.c$c.val, by=0.001)
exactCovProbVec <- Vectorize(exactCoverageProb, vectorize.args="c.val")
epsilon <- abs(exactCovProbVec(theta.diff=theta.diff, c.val=grid.c,
lambda=opt.lam.c$lambda, sigma.2=1, n=10)-0.95)
(new.c <- grid.c[which.min(epsilon)])


