REmod <-
function (X, y, ind, tind, n, k, t, nT, w, w2, coef0 = rep(0, 2),
    hess = FALSE, trace = trace, x.tol = 1.5e-18, rel.tol = 1e-15,
    ...)
{

## optimizing version 1:
    ##
    ## exploit ordering reversal
    ## and bdsmatrix functions as in ssrmod, sarsrmod, sarREmod...
    ##
    ## a) lag y etc.
    ## b) reverse ordering and exploit bds nature of vcov(e)
    ##
    ## maybe exploit analytical inverse of the submatrix block (gains on
    ## large-N problems)?? but how likely is it to have laaaarge T?

    ## extensive function rewriting, Giovanni Millo 04/10/2010
    ## structure:
    ## a) specific part
    ## - set names, bounds and initial values for parms
    ## - define building blocks for likelihood and GLS as functions of parms
    ## - define likelihood
    ## b) generic part(independent from ll.c() and #parms)
    ## - fetch covariance parms from max lik
    ## - calc last GLS step
    ## - fetch betas
    ## - calc final covariances
    ## - make list of results

    ## bdsmatrix and solve.bdsmatrix are now imported
    
    ## set names for final parms vectors
    nam.beta <- dimnames(X)[[2]]
    nam.errcomp <- c("phi")

    ## initialize values for optimizer
    myparms0 <- coef0
    ## set bounds for optimizer
    lower.bounds <- c(1e-08)
    upper.bounds <- c(1e+09)

    ## rearranging module
    ## save this for eventually re-rearranging output
    oo.0 <- order(tind, ind)
    ## reorder as stacked time series, as in std. panels
    oo <- order(ind, tind)
    X <- X[oo, ]
    y <- y[oo]
    #wy <- wy[oo]
    ind <- ind[oo]
    tind <- tind[oo]

    ## modules for likelihood
    bSigma <- function(phipsi, n, t, w) {
        ## single block of the original *scaled* covariance
        ## maintain w for homogeneity with generic part
        Jt <- matrix(1, ncol = t, nrow = t)
        It <- diag(1, t)
        ## retrieve parms
        phi <- phipsi[1]
        ## psi not used: here passing 2 parms, but works anyway
        ## because psi is last one
        ## calc inverse
        bSigma <- phi * Jt + It
        bSigma
    }
    detSigma <- function(phi, n, t) {
        detSigma <- -n/2 * log(t * phi + 1)
        detSigma
    }
    fullSigma <- function(phipsi, n, t, w) {
        sigma.i <- bSigma(phipsi, n, t, w)
        fullSigma <- bdsmatrix(rep(t, n), rep(as.numeric(sigma.i),
            n))
        fullSigma
    }

    ## likelihood function, both steps included
    ll.c <- function(phipsi, y, X, n, t, w, w2, wy) {
        ## retrieve parms
        phi <- phipsi[1]
        ## calc inverse sigma
        sigma <- fullSigma(phipsi, n, t, w)
        ## do GLS step to get e, s2e
        glsres <- GLSstepBDS(X, y, sigma)
        e <- glsres[["ehat"]]
        s2e <- glsres[["sigma2"]]
        ## calc ll
        due <- detSigma(phi, n, t)
        tre <- -n * t/2 * log(s2e)
        quattro <- -1/(2 * s2e) * crossprod(e, solve(sigma, e))
        const <- -(n * t)/2 * log(2 * pi)
        ll.c <- const + due + tre + quattro
        ## invert sign for minimization
        llc <- -ll.c
    }

    ## generic from here

    ## GLS step function for bdsmatrices
    GLSstepBDS <- function(X, y, sigma) {
        b.hat <- solve(crossprod(X, solve(sigma, X)), crossprod(X,
            solve(sigma, y)))
        ehat <- y - X %*% b.hat
        sigma2ehat <- crossprod(ehat, solve(sigma, ehat))/(n * t)
        return(list(betahat=b.hat, ehat=ehat, sigma2=sigma2ehat))
    }

    ## lag y unneeded here, keep parm for compatibility
    wy <- NULL

    ## max likelihood
    optimum <- nlminb(start = myparms0, objective = ll.c,
                      gradient = NULL, hessian = NULL,
                      y = y, X = X, n = n, t = t, w = w, w2 = w2, wy = wy,
                      scale = 1, control = list(x.tol = x.tol,
                                 rel.tol = rel.tol, trace = trace),
                      lower = lower.bounds, upper = upper.bounds)

    ## log likelihood at optimum (notice inverted sign)
    myll <- -optimum$objective
    ## retrieve optimal parms
    myparms <- optimum$par

    ## one last GLS step at optimal vcov parms
    sigma <- fullSigma(myparms, n, t, w)
    beta <- GLSstepBDS(X, y, sigma)

    ## final vcov(beta)
    covB <- as.numeric(beta[[3]]) *
        solve(crossprod(X, solve(sigma, X)))

    ## final vcov(errcomp)
    covTheta <- solve(-fdHess(myparms, function(x) -ll.c(x,
        y, X, n, t, w, w2, wy))$Hessian)          # lag-specific line: wy
    nvcovpms <- length(nam.errcomp)
    covAR <- NULL
    covPRL <- covTheta[1:nvcovpms, 1:nvcovpms, drop=FALSE]

    ## final parms
    betas <- as.vector(beta[[1]])
    sigma2 <- as.numeric(beta[["sigma2"]])    
    arcoef <- NULL
    errcomp <- myparms[which(nam.errcomp!="psi")]
    names(betas) <- nam.beta
    names(errcomp) <- nam.errcomp[which(nam.errcomp!="psi")]

    dimnames(covB) <- list(nam.beta, nam.beta)
    dimnames(covPRL) <- list(names(errcomp), names(errcomp))

    ## result
    RES <- list(betas = betas, arcoef=arcoef, errcomp = errcomp,
                covB = covB, covAR=covAR, covPRL = covPRL, ll = myll,
                sigma2 = sigma2)

    return(RES)
}
