# Method SOBI

SOBI <- function(X,...) UseMethod("SOBI")

# main function for SOBI
#
# input:
#  x = data matrix
#  k = either a intereger value, or a vector with integers. If a, then the lagged ACOVs 1:k are used otherwise those mentioned in k. 
#  method = method used for joint diagonalization, either rjd or djd
#  eps = eps for the JD
#  maxiter = maxiter for the JD

# output
#
# list of class "bss" with components
#   W = unmixing matrix
#   EV = Eigenvalues
#   k = lag used
#   S = sources as a time series object

SOBI.default <- function(X, k=12, method="frjd", eps = 1e-06, maxiter = 100, ...)
    {
    if (length(k)==1) k <- 1:k 
    nk <- length(k)
    method <- match.arg(method, c("rjd", "djd", "frjd"))
    
    MEAN <- colMeans(X)
    COV <- cov(X)
    EVD <- eigen(COV, symmetric = TRUE)
    COV.sqrt.i <- EVD$vectors %*% tcrossprod(diag(EVD$values^(-0.5)),EVD$vectors)
    X.C <- sweep(X,2,MEAN,"-")
    Y <- tcrossprod(X.C,COV.sqrt.i)
    p <- ncol(X)
    R <- array(0, dim=c(p,p,nk))
    n <- nrow(X) 
    
    for (i in 1:nk){
        Yt <- Y[1:(n-k[i]),]
        Yti <- Y[(1+k[i]):n,]
        Ri <- crossprod(Yt,Yti)/nrow(Yt)
        R[,,i] <- (Ri+t(Ri))/2
        }
    
    
    JD <- switch(method,
        frjd = {
                      frjd(R, eps = eps, maxiter = maxiter)$V
                      }
        ,
        "rjd"={
               rjd(R, eps=eps, maxiter=maxiter)$V
               }
        ,
        "djd"={
               djd(R, eps=eps, maxiter=maxiter,...)
               }
        )
    W <- crossprod(JD, COV.sqrt.i)
    #W <- diag(sign(rowMeans(W)))%*%W
    W <- sweep(W, 1, sign(rowMeans(W)), "*")
    S <- tcrossprod(X.C, W)
    
    acs <- acf(S, lag.max=max(k), plot=FALSE)$acf
    ssq_ac <- NULL
    for(j in 1:p){
      ssq_ac[j]<-sum(acs[k,j,j]^2)
    }
    ord <- order(ssq_ac, decreasing=TRUE)
    P <- matrix(0,p,p)
    for(j in 1:p){
      P[j,ord[j]]<-1
    }  
    S <- S[,ord]
    W <- P %*% W

    S <- ts(S, names=paste("Series",1:p))
    
    RES <- list(W=W, k=k, method=method, S=S)
    class(RES) <- "bss"
    RES
    }

SOBI.ts <- function(X, ...)
    {
    x <- as.matrix(X)
    RES <- SOBI.default(x,...)
    S <- RES$S
    attr(S, "tsp") <- attr(X, "tsp")
    RES$S <- S
    RES
    }
