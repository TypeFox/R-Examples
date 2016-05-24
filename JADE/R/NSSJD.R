# Method NSS.JD

NSS.JD <- function(X,...) UseMethod("NSS.JD")

# main function for NSS.JD
#
# input:
#  x = data matrix
#  K = number of intervals, default 12, obsolete when n.cuts is provided. 
#  Tau = lag of the autocovariance computed for each interval, default is 0 (= covariance matrix)
#  n.cuts = Cut of for intervals, must be of the form c(1, ..., n) where n is the sample size. If NULL, then K is used
#  eps = eps for the JD
#  maxiter = maxiter for the JD

# output
#
# list of class "bss" with components
#   W = unmixing matrix
#   k = lag used
#   n.cut = n.cut used
#   K = number of intervals used
#   S = sources as a time series object


NSS.JD.default <- function(X, K=12, Tau=0, n.cuts=NULL, eps = 1e-06, maxiter = 100, ...)
    {
    n <- nrow(X)
    MEAN <- colMeans(X)
    COV <- cov(X)
    EVD.COV <- eigen(COV, symmetric=TRUE)
    COV.sqrt.inv <- EVD.COV$vectors %*% tcrossprod(diag(sqrt(1/EVD.COV$values)), EVD.COV$vectors)
    X.C <- sweep(X,2,MEAN,"-")
    Y <- tcrossprod(X.C, COV.sqrt.inv)
    
    p <- ncol(X)
    
    if (is.null(n.cuts)) n.cuts <- ceiling(seq(1,n,length=K+1)) else K <- length(n.cuts)-1
    
    R <- array(0, dim=c(p,p,K))
    
    N.cuts <- n.cuts + c(rep(0,K),1)
    
    for (i in 1:K){
        R[,,i] <- M.x(Y[N.cuts[i]:(N.cuts[i+1]-1),], Tau=Tau)
        }
    
    W  <- crossprod(frjd(R, eps=eps, maxiter=maxiter)$V, COV.sqrt.inv) 
    S <- tcrossprod(X.C, W)
    S <- ts(S, names=paste("Series",1:p))
    
    RES <- list(W=W, k=Tau, n.cut=n.cuts, K=K, S=S)
    class(RES) <- "bss"
    RES
    }


NSS.JD.ts <- function(X, ...)
    {
    x <- as.matrix(X)
    RES <- NSS.JD.default(x,...)
    S <- RES$S
    attr(S, "tsp") <- attr(X, "tsp")
    RES$S <- S
    RES
    }
