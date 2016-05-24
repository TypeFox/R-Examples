shrink.estim <-
function(x,tar) {

    if (is.matrix(x)==TRUE && is.numeric(x)==FALSE) stop("The data matrix must be numeric!")

    p <- ncol(x) ; n <- nrow(x)
    covm <- cov(x) ; corm <- cor(x)
    xs <- scale(x,center=TRUE,scale=TRUE)
    
    v <- (n/((n-1)^3))*( crossprod(xs^2) - 1/n*(crossprod(xs))^2 )
    diag(v) <- 0
    
    m <- matrix(rep(apply(xs**2,2,mean),p),p,p)
    f <- (n/(2*(n-1)^3))*(  crossprod(xs**3,xs) + crossprod(xs,xs**3) - (m+t(m))*crossprod(xs)  )
    diag(f) <- 0 ; f[tar == 0] <- 0
    
    corapn <-  cov2cor(tar)
    d      <- (corm - corapn)^2
    lambda <- (sum(v)- sum(corapn*f))/sum(d)
    lambda <- max(min(lambda, 1), 0)
    shrink.cov <-lambda*tar+(1-lambda)*covm
    
    return(list(shrink.cov, c("The shrinkage intensity lambda is:",round(lambda,digits=4))))
}

