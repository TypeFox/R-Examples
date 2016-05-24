# Method NSS.JD

NSS.TD.JD <- function(X,...) UseMethod("NSS.TD.JD")

# main function for NSS.TD.JD
#
# input:
#  x = data matrix
#  K = number of intervals, default 12, obsolete when n.cuts is provided. 
#  Tau = lags used to compute the autocovariance matrices for each interval, default = 0:11
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


NSS.TD.JD.default <- function(X, K=12, Tau=0:11, n.cuts=NULL, eps = 1e-06, maxiter = 100, ...)
    {
    n <- nrow(X)
    MEAN <- colMeans(X)
    COV <- cov(X)
    EVD.COV <- eigen(COV, symmetric=TRUE)
    COV.sqrt.inv <- EVD.COV$vectors %*% tcrossprod(diag(sqrt(1/EVD.COV$values)),EVD.COV$vectors)
    X.C <-  sweep(X,2,MEAN,"-")
    Y <-  tcrossprod(X.C, COV.sqrt.inv)
    
    p <- ncol(X)
    
    if (is.null(n.cuts)) n.cuts <- ceiling(seq(1,n,length=K+1)) else K <- length(n.cuts)-1
    
    N.cuts <- n.cuts + c(rep(0,K),1)
    
    L<- length(Tau)
    
    R <- array(0, dim=c(p,p,L*K))
    
    ii<-1
    for (i in 1:K){
        Y.i<-Y[N.cuts[i]:(N.cuts[i+1]-1),]
            for (j in 1:L){
                R[,,ii] <- M.x(Y.i, Tau=Tau[j])
                ii<-ii+1
                }
        }
    
    W  <- crossprod(frjd(R, eps=eps, maxiter=maxiter)$V, COV.sqrt.inv)
    S <- tcrossprod(X.C,W)
    S <- ts(S, names=paste("Series",1:p))
    
    RES <- list(W=W, k=Tau, n.cut=n.cuts, K=K, S=S)
    class(RES) <- "bss"
    RES
    }


NSS.TD.JD.ts <- function(X, ...)
    {
    x <- as.matrix(X)
    RES <- NSS.TD.JD.default(x,...)
    S <- RES$S
    attr(S, "tsp") <- attr(X, "tsp")
    RES$S <- S
    RES
    }


##################################################################
# helper function M.x                                            #
# Autocovariance matrix for a centered time series at lag Tau    #
# (symmetrized)                                                  #
##################################################################

M.x <- function(X,Tau=0)
    {
    n<- nrow(X)
    Xt <- X[1:(n-Tau),]
    Xti <- X[(1+Tau):n,]
    Ri <- crossprod(Xt,Xti)/nrow(Xt)
    (Ri+t(Ri))/2 
    }
