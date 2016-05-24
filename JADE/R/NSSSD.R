# Method NSS.SD

NSS.SD <- function(X,...) UseMethod("NSS.SD")

# main function for NSS.SD
#
# input:
#  x = data matrix
#  k = a index at which point to cut the series, if NULL, then n.cut <- ceiling(n/2)


# output
#
# list of class "bss" with components
#   W = unmixing matrix
#   EV = Eigenvalues
#   n.cut = n.cut used
#   S = sources as a time series object


NSS.SD.default <- function(X, n.cut=NULL, ...)
    {
    n <- nrow(X)
    p <- ncol(X)
    MEAN <- colMeans(X)
    Y <- sweep(X,2,MEAN,"-") 
    
    if (!is.null(n.cut) & length(n.cut)==1)  n.cut <- c(1,n.cut,n)
    if (is.null(n.cut)) n.cut <- c(1,ceiling(n/2),n)
    if (length(n.cut)!=3) stop("'n.cut' must be NULL, or an integer or a vector of integers of length 3")
    Y1 <- Y[n.cut[1]:n.cut[2],]
    Y2 <- Y[(n.cut[2]+1):n.cut[3],]
    
    M.y1 <- cov(Y1)
    M.y2 <- cov(Y2)
    
    EVD.M.y1 <- eigen(M.y1, symmetric=TRUE)
    M.y1.sqrt.inv <- EVD.M.y1$vectors %*% tcrossprod(diag(sqrt(1/EVD.M.y1$values)), EVD.M.y1$vectors)
    
    M.12 <- M.y1.sqrt.inv %*% tcrossprod(M.y2, M.y1.sqrt.inv)
    
    EVD <- eigen(M.12,, symmetric=TRUE)
    
    W <- crossprod(EVD$vectors, M.y1.sqrt.inv)
    S <- tcrossprod(Y,W)
    S <- ts(S, names=paste("Series",1:p))
    
    RES <- list(W=W, EV=EVD$values, n.cut=n.cut, S=S)
    class(RES) <- "bss"
    RES    
    }



NSS.SD.ts <- function(X, ...)
    {
    x <- as.matrix(X)
    RES <- NSS.SD.default(x,...)
    S <- RES$S
    attr(S, "tsp") <- attr(X, "tsp")
    RES$S <- S
    RES
    }
