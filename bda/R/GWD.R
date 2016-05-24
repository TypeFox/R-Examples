pgwd <- function(q, alpha, sigma, lambda){
    stopifnot(sigma>0)
    x.neg <- q <= 0
    res <- rep(0, length(q))
    if(any(x.neg))
        res[x.neg] <- 0
    if(any(!x.neg)){
        q <- q[!x.neg]
        if(lambda == 0){
            u <- 1-exp(-(q/sigma)^(1/alpha))
            sele <- u >= 1
            if(any(sele)) u[sele] <- 1
            sele <- u < 0
            if(any(sele)) u[sele] <- 0
        }else{
            tmp <- 1 - lambda * (q/sigma)^(1/alpha)
            sele <- tmp <= 0
            if(any(sele)) tmp[sele] <- 0
            tmp <- 1-tmp^(1/lambda)
            tmp[sele] <- 1
        }
        res[!x.neg] <- tmp
    }
    res
}

dgwd <- function(x, alpha, sigma, lambda){
    stopifnot(sigma>0)
    res <- rep(0, length(x))
    x.finite <- is.finite(x)
    if(any(!x.finite))
        res[!x.finite] <- 0

    x.neg <- x <= 0
    if(any(x.neg))
        res[x.neg] <- 0

    sele <- !x.neg & x.finite
    if(any(sele)){
        q <- x[sele]
        u <- pgwd(q,alpha,sigma,lambda)
        if(lambda == 0){
            f <- (1-u)/alpha/qgwd(u,alpha-1,sigma,lambda)
        }else{
            f <- (1-u)^(1-lambda)/alpha/qgwd(u,alpha-1,sigma,lambda)
        }
        res[sele] <- f
    }
    res
}

qgwd <- function(p, alpha, sigma, lambda){
    stopifnot(p>=0 && p<=1)
    stopifnot(sigma>0)
    if(lambda == 0)
        res <- sigma*(-log(1-p))^alpha
    else
        res <- sigma*((1-(1-p)^lambda)/lambda)^alpha
    res
}

rgwd <- function(n, alpha, sigma, lambda){
    stopifnot(n>0)
    stopifnot(alpha>0&&sigma>0)
    stopifnot(is.finite(n))
    n <- round(n)
    u <- runif(n)
    qgwd(u, alpha=alpha, sigma=sigma,lambda=lambda)
}
