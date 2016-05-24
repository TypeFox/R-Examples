### R code from vignette source 'rvgtest.Rnw'

###################################################
### code chunk number 1: rvgtest.Rnw:151-155
###################################################
## load library
library("rvgtest")
## set seed (we want to have the same "random" plot in each session)
set.seed(123456789)


###################################################
### code chunk number 2: intro.Rnw:63-66 (eval = FALSE)
###################################################
## X <- rbeta(10000, shape1=3, shape2=0.05)
## U <- pbeta(X, shape1=3, shape2=0.05)
## hist(U)


###################################################
### code chunk number 3: hist-beta-1
###################################################
X <- rbeta(10000, shape1=3, shape2=0.05)
U <- pbeta(X, shape1=3, shape2=0.05)
hist(U)


###################################################
### code chunk number 4: intro.Rnw:86-89
###################################################
x <- c(1-2^(-52), 1-2^(-53), 1-2^(-54), 1-2^(-54))
pbeta(x, shape1=3, shape2=0.05)
1-x


###################################################
### code chunk number 5: histogram.Rnw:56-72
###################################################
##   Normal distribution
## .......................................................................
#
## Create frequency table.
## We use a sample of 20 times 10^5 random variates.
ft <- rvgt.ftable(n=1e5, rep=20, rdist=rnorm, qdist=qnorm)
#
## We can visualize this table by drawing the histogram
plot(ft)
#
## Perform a chi-square goodness-of-fit test
tr <- rvgt.chisq(ft)
tr

## We can visualize the p-values for increasing sample size
plot(tr)


###################################################
### code chunk number 6: hist-norm-1
###################################################
plot(ft)


###################################################
### code chunk number 7: pval-norm-1
###################################################
plot(tr)


###################################################
### code chunk number 8: histogram.Rnw:147-154
###################################################
##   Create and plot frequency table.
## .......................................................................
#
## Variant 1: use qnorm
ft1 <- rvgt.ftable(n=1e4, rep=10, rdist=rnorm, qdist=qnorm, plot=TRUE)
## Variant 3: use pnorm
ft2 <- rvgt.ftable(n=1e4, rep=10, rdist=rnorm, pdist=pnorm, plot=TRUE)


###################################################
### code chunk number 9: hist-norm-2a
###################################################
plot(ft1)


###################################################
### code chunk number 10: hist-norm-2b
###################################################
plot(ft2)


###################################################
### code chunk number 11: histogram.Rnw:187-199 (eval = FALSE)
###################################################
## ##   Perform goodness-of-fit tests
## ## .......................................................................
## #
## ## Create a frequency table.
## ft <- rvgt.ftable(n=1e5, rep=20, rdist=rnorm, qdist=qnorm)
## #
## ## Perform goodness-of-fit tests
## res.chisq <- rvgt.chisq(ft)
## res.mtest <- rvgt.Mtest(ft)
## #
## ## Show results in one plot
## plot.rvgt.htest(list(res.chisq, res.mtest))


###################################################
### code chunk number 12: pval-norm-2
###################################################
ft <- rvgt.ftable(n=1e5, rep=20, rdist=rnorm, qdist=qnorm)
res.chisq <- rvgt.chisq(ft)
res.mtest <- rvgt.Mtest(ft)
plot.rvgt.htest(list(res.chisq, res.mtest))


###################################################
### code chunk number 13: histogram.Rnw:225-243
###################################################
##   A buggy Gaussian random variate generator
## .......................................................................
#
## Use a buggy Gaussian random variate generator
RNGkind(normal.kind="Buggy Kinderman-Ramage")
#
## Create a frequency table.
ft <- rvgt.ftable(n=1e5, rep=50, rdist=rnorm, qdist=qnorm)
#
## Plot histogram (see Fig. 6)
plot(ft)
#
## Perform goodness-of-fit tests
res.chisq <- rvgt.chisq(ft)
res.mtest <- rvgt.Mtest(ft)
#
## Show results of tests (see Fig. 7)
plot.rvgt.htest(list(res.chisq, res.mtest))


###################################################
### code chunk number 14: histogram.Rnw:249-261 (eval = FALSE)
###################################################
## ##   Effect of increasing sample sizes
## ## .......................................................................
## #
## ## Create frequency table for buggy generator.
## RNGkind(normal.kind="Buggy Kinderman-Ramage")
## ft <- rvgt.ftable(n=1e5, rep=50, rdist=rnorm, qdist=qnorm)
## #
## ## Now plot histograms for the respective sample sizes 
## ## 1e5, 10*1e5, and 50*1e5.
## plot(ft,rows=1)
## plot(ft,rows=1:10)
## plot(ft,rows=1:50)


###################################################
### code chunk number 15: hist-buggynorm-1
###################################################
plot(ft)


###################################################
### code chunk number 16: pval-buggynorm-1
###################################################
res.chisq <- rvgt.chisq(ft)
res.mtest <- rvgt.Mtest(ft)
plot.rvgt.htest(list(res.chisq, res.mtest))


###################################################
### code chunk number 17: hist-buggynorm-2a
###################################################
plot(ft,rows=1)


###################################################
### code chunk number 18: hist-buggynorm-2b
###################################################
plot(ft,rows=1:10)


###################################################
### code chunk number 19: hist-buggynorm-2c
###################################################
plot(ft,rows=1:50)


###################################################
### code chunk number 20: inversion.Rnw:111-125 (eval = FALSE)
###################################################
## ##   u-error
## ## .......................................................................
## #
## ## Approximate inverse CDF of normal distribution using splines.
## aqn <- splinefun(x=pnorm((-100:100)*0.05), y=(-100:100)*0.05,
##                  method="monoH.FC")
## #
## ## Create a table of u-errors for the approximation errors.
## ## Use a sample of size of 10^5 random variates and 
## ## a resolution of only 100 points.
## ue <- uerror(n=1e5, res=100, aqdist=aqn, pdist=pnorm)
## #
## ## Plot u-errors.
## plot(ue)


###################################################
### code chunk number 21: uerror-1
###################################################
aqn <- splinefun(x=pnorm((-100:100)*0.05), y=(-100:100)*0.05,
                 method="monoH.FC")
ue <- uerror(n=1e5, res=100, aqdist=aqn, pdist=pnorm)
plot(ue)


###################################################
### code chunk number 22: uerror-2
###################################################
aqn <- splinefun(x=pnorm((-100:100)*0.05), y=(-100:100)*0.05,
                 method="monoH.FC")
ue <- uerror(n=1e5, res=1000, aqdist=aqn, pdist=pnorm)
plot(ue)


###################################################
### code chunk number 23: inversion.Rnw:159-171 (eval = FALSE)
###################################################
## ##   u-error
## ## .......................................................................
## #
## ## Approximate inverse CDF of normal distribution using splines.
## aqn <- splinefun(x=pnorm((-100:100)*0.05), y=(-100:100)*0.05,
##                  method="monoH.FC")
## #
## ## Create a table of u-errors for the approximation errors 
## ## for the subdomain (0.6, 0.65).
## ## Use a sample of size of 10^5 random variates.
## ue <- uerror(n=1e5, aqdist=aqn, pdist=pnorm, udomain=c(0.6,0.65))
## plot(ue)


###################################################
### code chunk number 24: uerror-3
###################################################
aqn <- splinefun(x=pnorm((-100:100)*0.05), y=(-100:100)*0.05,
                 method="monoH.FC")
ue <- uerror(n=1e5, aqdist=aqn, pdist=pnorm, udomain=c(0.6,0.65), plot=TRUE)


###################################################
### code chunk number 25: inversion.Rnw:191-206 (eval = FALSE)
###################################################
## ##   u-error
## ## .............................................................
## #
## ## Approximate inverse CDF of normal and gamma distribution using splines.
## aqn <- splinefun(x=pnorm((-100:100)*0.05), y=(-100:100)*0.05,
##                  method="monoH.FC")
## aqg <- splinefun(x=pgamma((0:500)*0.1,shape=5),
##                  y=(0:500)*0.1, method="monoH.FC")
## #
## ## Compute u-errors for these approximations
## uen <- uerror(n=1e5, aqdist=aqn, pdist=pnorm)
## ueg <- uerror(n=1e5, aqdist=aqg, pdist=pgamma, shape=5)
## #
## ## Plot u-errors
## plot.rvgt.ierror(list(uen,ueg), tol=1.e-6)


###################################################
### code chunk number 26: uerror-4
###################################################
aqn <- splinefun(x=pnorm((-100:100)*0.05), y=(-100:100)*0.05,
                 method="monoH.FC")
uen <- uerror(n=1e5, aqdist=aqn, pdist=pnorm)
aqg <- splinefun(x=pgamma((0:500)*0.1,shape=5),
                 y=(0:500)*0.1, method="monoH.FC")
ueg <- uerror(n=1e5, aqdist=aqg, pdist=pgamma, shape=5)
plot.rvgt.ierror(list(uen,ueg), tol=1.e-6)


###################################################
### code chunk number 27: inversion.Rnw:235-247 (eval = FALSE)
###################################################
## ##   Absolute x-error
## ## .............................................................
## #
## ## Approximate inverse CDF of normal distribution using splines.
## aqn <- splinefun(x=pnorm((-100:100)*0.05), y=(-100:100)*0.05,
##                  method="monoH.FC")
## #
## ## Create a table of absolute x-errors.
## xea <- xerror(n=1e5, aqdist=aqn, qdist=qnorm, kind="abs")
## #
## ## Plot x-errors
## plot(xea)


###################################################
### code chunk number 28: xerror-abs
###################################################
aqn <- splinefun(x=pnorm((-100:100)*0.05), y=(-100:100)*0.05, method="monoH.FC")
xea <- xerror(n=1e5, aqdist=aqn, qdist=qnorm, kind="abs")
plot(xea)


###################################################
### code chunk number 29: inversion.Rnw:263-275 (eval = FALSE)
###################################################
## ##   Relative x-error
## ## .......................................................................
## #
## ## Approximate inverse CDF of normal distribution using splines.
## aqn <- splinefun(x=pnorm((-100:100)*0.05), y=(-100:100)*0.05,
##                  method="monoH.FC")
## #
## ## Create a table of relative x-errors.
## xer <- xerror(n=1e5, aqdist=aqn, qdist=qnorm, kind="rel")
## #
## ## Plot x-errors
## plot(xer)


###################################################
### code chunk number 30: xerror-rel
###################################################
aqn <- splinefun(x=pnorm((-100:100)*0.05), y=(-100:100)*0.05, method="monoH.FC")
xer <- xerror(n=1e5, aqdist=aqn, qdist=qnorm, kind="rel")
plot(xer)


