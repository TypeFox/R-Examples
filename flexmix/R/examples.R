#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: examples.R 4859 2012-12-18 08:42:33Z gruen $
#

ExNPreg = function(n)
{
    if(n %% 2 != 0) stop("n must be even")
    
    x <- runif(2*n, 0, 10)
    mp <- exp(c(2-0.2*x[1:n], 1+0.1*x[(n+1):(2*n)]))
    mb <- binomial()$linkinv(c(x[1:n]-5, 5-x[(n+1):(2*n)]))

    data.frame(x=x,
               yn=c(5*x[1:n], 40-(x[(n+1):(2*n)]-5)^2)+3*rnorm(n),
               yp=rpois(2*n, mp),
               yb=rbinom(2*n, size=1, prob=mb),
               class = rep(1:2, c(n,n)),
               id1 = factor(rep(1:n, rep(2, n))),
               id2 = factor(rep(1:(n/2), rep(4, n/2))))
}


    
ExNclus = function(n=100)
{
    if(n %% 2 != 0) stop("n must be even")

    rbind(mvtnorm::rmvnorm(n, mean=rep(0,2)),
          mvtnorm::rmvnorm(n, mean=c(8,0), sigma=diag(1:2)),
          mvtnorm::rmvnorm(1.5*n, mean=c(-2,6), sigma=diag(2:1)),
          mvtnorm::rmvnorm(2*n, mean=c(4,4), sigma=matrix(c(1,.9,.9,1), 2)))
}

    
ExLinear <- function(beta, n, xdist="runif", xdist.args=NULL,
                     family=c("gaussian", "poisson"), sd=1, ...)
{
    family <- match.arg(family)
    
    X <- NULL
    y <- NULL
    k <- ncol(beta)
    d <- nrow(beta)-1

    n <- rep(n, length.out=k)
    if(family=="gaussian") sd <- rep(sd, length.out=k)
    xdist <- rep(xdist, length.out=d)
    
    if(is.null(xdist.args)){
        xdist.args <- list(list(...))
    }
    if(!is.list(xdist.args[[1]]))
        xdist.args <- list(xdist.args)
    
    xdist.args <- rep(xdist.args, length.out=d)
    
    for(i in 1:k)
    {
        X1 <- 1
        for(j in 1:d){
            xdist.args[[j]]$n <- n[i]
            X1 <- cbind(X1, do.call(xdist[j], xdist.args[[j]]))
        }

        X <- rbind(X, X1)
        xb <- X1 %*% beta[,i,drop=FALSE]
        if(family=="gaussian")
            y1 <- xb + rnorm(n[i], sd=sd[i])
        else
            y1 <- rpois(n[i], exp(xb))
    
        y <- c(y, y1)
    }
    X <- X[,-1,drop=FALSE]
    colnames(X) <- paste("x", 1:d, sep="")
               
    z <- data.frame(y=y, X)
    attr(z, "clusters") <- rep(1:k, n)
    z
}

