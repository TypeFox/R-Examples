CovMest <- function(x, r = 0.45, arp = 0.05, eps=1e-3, maxiter=120, control, t0, S0, initcontrol)
{

    ## Analize and validate the input parameters ...

    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.
    if(!missing(control)){
        defcontrol <- new("CovControlMest")      # default control
        if(r == defcontrol@r)       r <- control@r
        if(arp == defcontrol@arp)   arp <- control@arp
        if(eps == defcontrol@eps)   eps <- control@eps
        if(maxiter == defcontrol@maxiter) maxiter <- control@maxiter
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

    call <- match.call()
    method <- "M-Estimates"

    ## If not provided initial estimates, compute them as MVE
    ##  Take the raw estimates and standardise the covariance
    ##  matrix to determinant=1
    if(missing(t0) || missing(S0)){
        if(missing(initcontrol))
            init <- CovMve(x)
        else
            init <- restimate(initcontrol, x)
        cl <- class(init)
        if(cl == "CovMve" || cl == "CovMcd" || cl == "CovOgk"){
            t0 <- init@raw.center
            S0 <- init@raw.cov
            detS0 <-det(S0)
            detS02 <- detS0^(1.0/p)
            S0 <- S0/detS02
        } else{
            t0 <- init@center
            S0 <- init@cov
        }
    }

    ## calculate the constants M and c
    ## for the translated biweight function
    psix <- new("PsiBwt", n=n, p=p, r=r, alpha=arp)
    psix <- csolve(psix)
    mest <- iterM(psix, x, t0, S0, eps=eps, maxiter=maxiter)

    mah <- mahalanobis(x, mest$t1, mest$s)
    crit <- determinant(mest$s, logarithm = FALSE)$modulus[1]

    new("CovMest", call=call, cov=mest$s, center=mest$t1, n.obs=n,
        mah=mah, method=method, X=x, iter=mest$iter, crit=crit,
        wt=mest$wt, vt=mest$vt)
}
