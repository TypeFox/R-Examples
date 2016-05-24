##' @export
pbvn <- function(upper,rho,sigma) {
    if (!missing(sigma)) {
        rho <- cov2cor(sigma)[1,2]
        upper <- upper/diag(sigma)^0.5
    }
    arglist <- list("bvncdf",
                    a=upper[1],
                    b=upper[2],
                    r=rho,
                    PACKAGE="mets")
    res <- do.call(".Call",arglist)
    return(res)
}


## lower <- rbind(c(0,-Inf),(-Inf,0))
## upper <- rbind(c(Inf,0),(0,Inf))
## mu <- rbind(c(1,1),c(-1,1))
## sigma <- diag(2)+1
## pmvn(lower=lower,upper=upper,mu=mu,sigma=sigma)

##' @export
pmvn <- function(lower,upper,mu,sigma,cor=FALSE) {
    if (missing(sigma)) stop("Specify variance matrix 'sigma'")
    if (missing(lower)) {
        if (missing(upper)) stop("Lower or upper integration bounds needed")
        lower <- upper; lower[] <- -Inf
    }
    p <- ncol(rbind(lower))
    if (missing(upper)) {
        upper <- lower; upper[] <- Inf
    }
    if (missing(mu)) mu <- rep(0,p)
    sigma <- rbind(sigma)
    ncor <- p*(p-1)/2
    if (ncol(sigma)!=p && ncol(sigma)!=ncor)
        stop("Incompatible dimensions of mean and variance")    
    if (ncol(rbind(lower))!=p || ncol(rbind(upper))!=p)
        stop("Incompatible integration bounds")    
    arglist <- list("pmvn",
                    lower=rbind(lower),
                    upper=rbind(upper),
                    mu=rbind(mu),
                    sigma=rbind(sigma),
                    cor=as.logical(cor[1]))
    res <- do.call(".Call",arglist)
    return(as.vector(res))
}
