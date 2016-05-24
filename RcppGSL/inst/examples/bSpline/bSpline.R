
## This example illustrated use of RcppGSL using the 'Rcpp attributes' feature
##
## The example comes from Section 39.7 of the GSL Reference manual, and constructs
## a data set from the curve y(x) = \cos(x) \exp(-x/10) on the interval [0, 15] with
## added Gaussian noise --- which is then fit via linear least squares using a cubic
## B-spline basis functions with uniform breakpoints.
##
## Obviously all this could be done in R too as R can both generate data, and fit
## models including (B-)splines. But the point to be made here is that we can very
## easily translate a given GSL program (thanks to RcppGSL), and get it into R with
## ease thanks to Rcpp and Rcpp attributes.

require(Rcpp)                           # load Rcpp

sourceCpp("bSpline.cpp")                # compile two functions

dat <- genData()                        # generate the data
fit <- fitData(dat)                     # fit the model, returns matrix and gof measures

X <- fit[["X"]]                         # extract vectors
Y <- fit[["Y"]]

op <- par(mar=c(3,3,1,1))
plot(dat[,"x"], dat[,"y"], pch=19, col="#00000044")
lines(X, Y, col="orange", lwd=2)
par(op)
