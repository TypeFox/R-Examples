WRegEst <- function(x, y, delta, LS=TRUE, tau=0.5) {
# This R function computes the case weighted regression 
# estimator for randomly right censored data. It can compute the
# least squares estimator or quantile regression estimator. In the
# later case, it calls a function in the quantreg() package.
# Input:
# x is a matrix of N rows (covariates).
# y is the observed (censored) responses --- a vector of length N.
# delta is a vector of length N. delta =1 means (y) is not censored.
#           delta = 0 means y is right censored, i.e. the true
#        response is larger than y.
# LS  indicates if this is a least square regression or quantile regression
# tau if LS=TRUE then this is ignored, otherwise tau is used in rqfit
#
# Output:
# the estimate, \hat beta. 

n <- length(y)
x <- as.matrix(x)
xdim <- dim(x)
if ( xdim[1] != n ) stop("check dim of x")
if ( length(delta) != n ) stop("check length of delta")
if(any((delta!=0)&(delta!=1)))
   stop("delta must be 0(right-censored) or 1(uncensored)")


temp <- WKM(x=y, d=delta, zc=1:n)
KMweight <- temp$jump
norder <- order(y, -delta)
### the zero weights should be removed
Wplace <- which(KMweight > 0)
KMweight <- KMweight[Wplace]
ZZ <- y[norder]
ZZ <- ZZ[Wplace]
XX <- as.matrix(x[norder,])
xmat <- as.matrix(XX[Wplace,])
if(LS) estim <- coef(lm.wfit(x=xmat, y=ZZ, w=KMweight))
if(!LS) estim <- coef(rq.wfit(x=xmat, y=ZZ, weights=KMweight, tau=tau))
return(estim)
}
