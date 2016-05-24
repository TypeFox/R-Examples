covMest <- function(x, cor=FALSE,  r = 0.45, arp = 0.05, eps=1e-3, maxiter=120, control, t0, S0)
{
    .Deprecated(new="CovMest")

    ## Analize and validate the input parameters ...

    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.
    if(!missing(control)){
        defcontrol <- rrcov.control()      # default control
#        if(r == defcontrol$r)       r <- control$r
#        if(arp == defcontrol$arp)   arp <- control$arp
#        if(eps == defcontrol$eps)   eps <- control$eps
#        if(maxiter == defcontrol$maxiter) maxiter <- control$maxiter
    }

    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
            dimnames = list(names(x), deparse(substitute(x))))

    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok, , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
        stop("All observations have missing values!")
    dimn <- dimnames(x)
    n <- dx[1]
    p <- dx[2]
    if(n < 2 * p)
        stop("Need at least 2*(number of variables) observations ")

    ans <- list(method = "M-Estimates", call = match.call())

    ## If not provided initial estimates, compute them as MVE
    ##  Take the raw estimates and standardise the covariance
    ##  matrix to determinant=1
    if(missing(t0) || missing(S0)){
        mcd <- CovMve(x)
        t0 <- mcd@raw.center
        S0 <- mcd@raw.cov
        detS0 <-det(S0)
        detS02 <- detS0^(1.0/p)
        S0 <- S0/detS02

    }

    ## calculate the constants M and c
    ## for the translated biweight function
    psix <- new("PsiBwt", n=n, p=p, r=r, alpha=arp)
    psix <- csolve(psix)
    mest <- iterM(psix, x, t0, S0, eps=1e-3, maxiter=maxiter)

    ## this was the version without OO
    ##const <- csolve.bt(n, p, r, arp)
    ##mest <- .iterM(x, t0, S0, const$c1, const$M, eps, maxiter)

    ans$n.obs <- n
    ##ans$c1 <- const$c1
    ##ans$M <- const$M

    ans$c1 <- psix@c1
    ans$M <- psix@M

    ans$iter <- mest$iter
    ans$cov <- mest$s
    ans$center <- mest$t1
    ans$mah <- mahalanobis(x, mest$t1, mest$s)
    ans$crit <- determinant(mest$s, logarithm = TRUE)$modulus[1]
    if(cor && !is.null(ans$cov))
        cor <- cov2cor(ans$cov)

    class(ans) <- c("mest", "mcd")
    attr(ans, "call") <- sys.call()
    ans$method <- paste("M-Estimator.")
    ans$X <- x

    return(ans)
}
