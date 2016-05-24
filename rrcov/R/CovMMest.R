CovMMest <- function(x,
                    bdp=0.5,
                    eff=0.95,
                    eff.shape=TRUE,
                    maxiter=50,
                    trace=FALSE,
                    tolSolve=1e-7,
                    control
                )
{
    ## NOTES:
    ##  - in the functions rho, psi, and scaledpsi=psi/u (i.e. the weight function)
    ##      is used |x| <= c1
    ##
    ##  The bisquare rho function:
    ##
    ##              |   x^2/2 - x^4/2*c1^2 + x^6/6*c1^4         |x| <=  c1
    ##  rho(x) =    |
    ##              |   c1^2/6                                  |x| > c1
    ##
    rho <- function(u, cc)
    {
        w <- abs(u) <= cc
        v <- (u^2/2 * (1 - u^2/cc^2 + u^4/(3*cc^4))) * w + (1-w) * (cc^2/6)
        v
    }

    ##  The corresponding psi function: psi = rho'
    ##
    ##              |   x - 2x^3/c1^2 + x^5/c1^4         |x| <=  c1
    ##  psi(x) =    |
    ##              |   0                                |x| > c1
    ##
    ##      using ifelse is 3 times slower
    psi <- function(u, c1)
    {
        ##ifelse(abs(u) < c1, u - 2 * u^3/c1^2 + u^5/c1^4, 0)
        pp <- u - 2 * u^3/c1^2 + u^5/c1^4
        pp*(abs(u) <= c1)
    }

    ## weight function = psi(u)/u
    scaledpsi <- function(u, cc)
    {
        ##ifelse(abs(xx) < c1, xx - 2 * xx^3/c1^2 + xx^5/c1^4, 0)
        pp <- (1 - (u/cc)^2)^2
        pp <- pp * cc^2/6
        pp*(abs(u) <= cc)
    }

    ## the objective function, we solve loss.S(u, s, cc) = b for "s"
    loss <- function(u, s, cc)
        mean(rho(u/s, cc))

    ## Returns square root of the mahalanobis distances of x with respect to mu and sigma
    ## Seems to be somewhat more efficient than sqrt(mahalanobis()) - by factor 1.4!
    resdis <- function(x, mu, sigma)
    {
        central <- t(x) - mu
        sqdis <- colSums(solve(sigma, central) * central)
        dis <- sqdis^(0.5)
        dis
    }

    ##
    ## Compute weighted mean and covariance matrix
    ##  The covariance is scaled to have det=1
    ##
    covw <- function(x, wt=rep(1, nrow(x)))
    {
        if (is.data.frame(x))
            x <- as.matrix(x)
        else if (!is.matrix(x))
            stop("'x' must be a matrix or a data frame")
        if (!all(is.finite(x)))
            stop("'x' must contain finite values only")

        n <- nrow(x)
        p <- ncol(x)
        if(with.wt <- !missing(wt))
        {
            if(length(wt) != n)
                stop("length of 'wt' must equal the number of rows in 'x'")
            if(any(wt < 0) || (s <- sum(wt)) == 0)
                stop("weights must be non-negative and not all zero")
        }

        center <- colSums(wt * x)/s
        x <- sqrt(wt) * sweep(x, 2, center, check.margin = FALSE)
        cov <- crossprod(x)
        cov <- det(cov)^(-1/p) * cov
        ret <- list(cov = cov, center = center, n.obs = n)
        ret
    }


    ## Analize and validate the input parameters ...
    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.
    scontrol <- CovControlSest()
    if(!missing(control)){
        defcontrol <- CovControlMMest()     # default control
        if(bdp == defcontrol@bdp)           bdp <- control@bdp          #
        if(eff == defcontrol@eff)           eff <- control@eff          #
        if(maxiter == defcontrol@maxiter)   maxiter <- control@maxiter  # for bisquare and rocke

        if(trace == defcontrol@trace)       trace <- control@trace
        if(tolSolve == defcontrol@tolSolve) tolSolve <- control@tolSolve
        scontrol <- control@sest
        scontrol@bdp = bdp
    }

    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
                    dimnames = list(names(x), deparse(substitute(x))))

    xcall <- match.call()

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

    ## compute the constant c1 for the rho function of the MM-estimates
    c1 <- .csolve.bw.MM(p, eff, eff.shape=eff.shape)

    if(trace)
        cat("\nMM-EST...: bdp, eff, eff.shape, p, c1=", bdp, eff, eff.shape, p, c1, "\n")

    ## Compute the initial S-estimate
    ss <- CovSest(x, control=scontrol, trace=trace)
    if(trace)
    {
        cat("\nInitial S-estimates. Scale=", ss@crit, "\n")
        print(ss)
    }

    scale <- ss@crit
    cova <- getShape(ss)
    center <- getCenter(ss)

    rdis <- resdis(x, center, cova)
    keep <- ocrit <- crit <- loss(rdis, scale, c1)
    if(trace)
    {
        cat("\nInitial robust distances:\n")
        print(rdis)
    }

    ## Compute M-estimate with auxiliary scale.
    ##  Start from an S-estimate and perform reweighted
    ##  location/covariance steps
    for(iter in 1:maxiter)
    {
        if(trace)
            cat("\nITER=", iter,"\n")

        w <- scaledpsi(rdis/scale, c1)

        center <- as.vector(crossprod(w, x) / sum(w))
        ccc <- covw(x, w)
        center <- ccc$center
        cova <- ccc$cov
        rdis <- resdis(x, center, cova)
        crit <- loss(rdis, scale, c1)

##        cat("\noldobj, newobj, iter:",oldobj, newobj, oldobj-newobj, oldobj - newobj > tolSolve, iter,"\n")
        if(ocrit - crit <= tolSolve)
            break
        ocrit <- crit
    }
    cova <- cova*scale^2

    ## Check if the loss function was at all reduced by
    ##  the M-iterations. If not (which should never happen)
    ##  simply return the S-estimate
    if(crit > keep)
    {
        center <- ss@center
        cova <- ss@cov
    }

    ans <- new("CovMMest",
               call = xcall,
               iter=iter,
               crit=ss@crit,
               cov=cova,
               center=center,
               c1=c1,
               n.obs=n,
               X = as.matrix(x),
               sest=ss,
               method="MM-estimates")
    ans
}

##
## Compute the constant for the second Tukey Biweight rho-function for MM
##  with for fixed shape-efficiency
##
## Adapted from Gert Willems:
##      http://users.ugent.be/~svaelst/software/MMPCAboot.html
##
.csolve.bw.MM <- function(p, eff, eff.shape=TRUE)
{
    ## (taken from Claudia Becker's Sstart0 program)
    chi.int <- function(p, a, c1)
    {
        ## partial expectation d in (0,c1) of d^a under chi-squared p
        return(exp(lgamma((p + a) /2) - lgamma(p/2)) * 2^{a/2} * pchisq(c1^2, p + a))
    }

    loceff.bw <- function(p, c1)
    {
        # called by csolve.bw.MM(); computes location efficiency corresponding to c1
        alpha1 <- 1/p * (chi.int(p,2,c1) - 4*chi.int(p,4,c1)/(c1^2) + 6*chi.int(p,6,c1)/(c1^4) - 4*chi.int(p,8,c1)/(c1^6) + chi.int(p,10,c1)/(c1^8))
        beta1.1 <- chi.int(p,0,c1) - 2*chi.int(p,2,c1)/(c1^2) + chi.int(p,4,c1)/(c1^4)
        beta1.2 <- chi.int(p,0,c1) - 6*chi.int(p,2,c1)/(c1^2) + 5*chi.int(p,4,c1)/(c1^4)
        beta1 <- (1-1/p)*beta1.1 + 1/p*beta1.2
        return( beta1^2 / alpha1 )
    }


    sigma1.bw <- function(p, c1)
    {
        gamma1.1 <- chi.int(p,2,c1) - 6*chi.int(p,4,c1)/(c1^2) + 5*chi.int(p,6,c1)/(c1^4)
        gamma1.2 <- chi.int(p,2,c1) - 2*chi.int(p,4,c1)/(c1^2) + chi.int(p,6,c1)/(c1^4)
        gamma1 <- (gamma1.1 + (p+1)*gamma1.2) / (p+2)
        sigma1.0 <- chi.int(p,4,c1) - 4*chi.int(p,6,c1)/(c1^2) + 6*chi.int(p,8,c1)/(c1^4) - 4*chi.int(p,10,c1)/(c1^6) + chi.int(p,12,c1)/(c1^8)
        return(sigma1.0 / (gamma1^2) * p/(p+2))
    }

    maxit <- 1000
    eps <- 1e-8

    ## ctest <- csolve.bw.asymp(p,.5)
    cold <- ctest <- -.4024 + 2.2539 * sqrt(p) # very precise approximation for c corresponding to 50% bdp
    for(iter in 1:maxit)
    {
        ctest <- if(eff.shape) ctest * eff * sigma1.bw(p, ctest) else ctest * eff / loceff.bw(p, ctest)
        if(abs(cold-ctest) <= eps)
            break
        cold <- ctest
    }
    return(ctest)
}
