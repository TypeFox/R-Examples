# Method AMUSE

AMUSE <- function(x,...) UseMethod("AMUSE")

# main function for AMUSE
#
# input:
#  x = data matrix
#  k = lag, default=1
#
# output
#
# list of class "bss" with components
#   W = unmixing matrix
#   EV = Eigenvalues
#   k = lag used
#   S = sources as a time series object

AMUSE.default <- function(x,k=1,...) 
    {
    na.fail(x)
    p <- ncol(x)
    n <- nrow(x)
    
    x.c <- sweep(x,2,colMeans(x), "-")
    COV <- crossprod(x.c)/(n-1)
    COV.EVD <- eigen(COV, symmetric = TRUE)
    COV.sqrt.i <- COV.EVD$vectors %*% tcrossprod(diag(COV.EVD$values^(-0.5)), COV.EVD$vectors)

    ACOVk <- crossprod(x.c[1:(n-k),],x.c[(k+1):n,])/(n-k)
    ACOVk.sym <- (ACOVk + t(ACOVk))/2
    COViACOVk <- crossprod(COV.sqrt.i,tcrossprod(ACOVk.sym,COV.sqrt.i))
    
    EVD <- eigen(COViACOVk, symmetric=TRUE)
    
    W <- crossprod(EVD$vectors, COV.sqrt.i)
    W <- sweep(W, 1, sign(rowMeans(W)), "*")

    EV <- EVD$values
    S <- tcrossprod(x,W)
    S <- ts(S, names=paste("Series",1:p))
    
    RES <- list(W=W, EV=EV, k=k, S=S)
    class(RES) <- "bss"
    RES
    }
    
# Function for AMUSE when input is a multivariate time series

AMUSE.ts<- function(x,...)
    {
    X <- as.matrix(x)
    RES <- AMUSE.default(X,...)
    S <- RES$S
    attr(S, "tsp") <- attr(x, "tsp")
    RES$S <- S
    RES
    }

# printing method for objects of class "bss"

`print.bss` <-
function(x, ...)
    {
    print.listof(x[names(x)!="S"], ...)
    }

# ploting method for objects of class "bss"
# is R's basic time series plot

`plot.bss` <-
function(x, ...)
    {
    S <- x$S
    
    if (is.ts(S)) plot.ts(S,y=NULL, ...) else pairs(S, ...)
    }

# coef method for objects of class "bss"
# extracts the unmixing matrix

`coef.bss` <-
function(object, ...)
    {
    object$W
    }

# extracting the estimated signals

`bss.components` <-
function(object)
    {
    object$S
    }
