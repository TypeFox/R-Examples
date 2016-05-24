#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: distances.R 3 2013-06-12 10:06:43Z leisch $
#

distEuclidean <- function(x, centers)
{
    if(ncol(x)!=ncol(centers))
        stop(sQuote("x")," and ",sQuote("centers"),
             " must have the same number of columns")
    z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
    for(k in 1:nrow(centers)){
        z[,k] <- sqrt( colSums((t(x) - centers[k,])^2) )
    }
    z
}

distManhattan <- function(x, centers)
{
    if(ncol(x)!=ncol(centers))
        stop(sQuote("x")," and ",sQuote("centers"),
             " must have the same number of columns")
    z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
    for(k in 1:nrow(centers)){
        z[,k] <- colSums(abs(t(x) - centers[k,]))
    }
    z
}

distMax <- function(x, centers)
{
    if(ncol(x)!=ncol(centers))
        stop(sQuote("x")," and ",sQuote("centers"),
             " must have the same number of columns")
    z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
    for(k in 1:nrow(centers)){
        z[,k] <- apply(abs(t(x) - centers[k,]), 2, max)
    }
    z
}

distJaccard <- function(x, centers)
{
    if(ncol(x)!=ncol(centers))
        stop(sQuote("x")," and ",sQuote("centers"),
             " must have the same number of columns")
    xc <- x %*% t(centers)
    nenner <-
        matrix(rowSums(x), nrow=nrow(x), ncol=nrow(centers)) +
            matrix(rowSums(centers), nrow=nrow(x), ncol=nrow(centers),
                   byrow=TRUE) - xc

    z <- 1 - xc/nenner
    z[nenner<sqrt(.Machine$double.eps)] <- 0
    z
}

distCanberra <- function(x, centers)
{
    if(ncol(x)!=ncol(centers))
        stop(sQuote("x")," and ",sQuote("centers"),
             " must have the same number of columns")
    z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
    tx <- t(x)
    for(k in 1:nrow(centers)){
        d <- abs(tx-centers[k,])
        s <- abs(tx+centers[k,])
        q <- d/s
        q[s<.Machine$double.eps] <- 0
        ## in dist() erhöhen doppelte nullen die distanz um einen
        ## faktor -> abgekupfert für konsistenz. 
        z[,k] <- colSums(q) * ncol(x) / colSums(s>.Machine$double.eps)
    }
    z
}

distMinkowski <- function(x, centers, p=2)
{
    if(ncol(x)!=ncol(centers))
        stop(sQuote("x")," and ",sQuote("centers"),
             " must have the same number of columns")
    z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
    for(k in 1:nrow(centers)){
        z[,k] <- colSums(abs(t(x) - centers[k,])^p)^(1/p)
    }
    z
}

dist2 <- function(x, y, method = "euclidean", p=2){

    if(any(is.na(x)) || any(is.na(y)))
        stop("Cannot handle missing values!")
    
    x <- as(x, "matrix")

    if(is.vector(y) && (length(y)<=ncol(x)))
        y <- matrix(y, nrow=1, ncol=ncol(x))
    else
        y <- as(y, "matrix")
    
    METHODS <- c("euclidean", "maximum", "manhattan", "canberra", 
                 "binary", "minkowski")
    method <- match.arg(method, METHODS)

    z <- switch(method,
                euclidean = distEuclidean(x, y),
                maximum = distMax(x, y),
                manhattan = distManhattan(x, y),
                canberra = distCanberra(x, y),
                binary = distJaccard(x!=0, y!=0),
                minkowski = distMinkowski(x, y, p=p))
    rownames(z) <- rownames(x)
    colnames(z) <- rownames(y)
    z
}
           



###**********************************************************

centMean <- function(x) colMeans(x)

centMedian <- function(x) apply(x, 2, median)

centOptim <- function(x, dist)
{
    foo <- function(p)
        sum(dist(x, matrix(p, nrow=1)))

    optim(colMeans(x), foo)$par
}

centOptim01 <- function(x, dist)
{
    foo <- function(p)
        sum(dist(x, matrix(p, nrow=1)))

    optim(colMeans(x), foo, lower=0, upper=1, method="L-BFGS-B")$par
}

###**********************************************************

distAngle <- function(x, centers)
{
    if(ncol(x)!=ncol(centers))
        stop(sQuote("x")," and ",sQuote("centers"),
             " must have the same number of columns")
    z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
    for(k in 1:nrow(centers)){
        if(any(is.na(centers[k,])))
            z[,k] <- - Inf
        else
            z[,k] <- apply(x, 1, function(a) sum(a*centers[k,]))
    }
    1-z
}


centAngle <- function(x)
{
    z <- colMeans(x)
    z/sqrt(sum(z^2))
}

wcentAngle <- function(x, weights)
{
    z <- colMeans(x*normWeights(weights))
    z/sqrt(sum(z^2))
}

###**********************************************************

distCor <- function(x, centers)
{
   z <- matrix(0,nrow(x), ncol=nrow(centers))
   for(k in 1:nrow(centers)){
      z[,k] <- 1 - cor(t(x), centers[k,])
   }
   z
}   

