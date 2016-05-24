ssrREmod <-
function (X, y, ind, tind, n, k, t, nT, w, w2, coef0 = rep(0, 3),
    hess = FALSE, trace = trace, x.tol = 1.5e-18, rel.tol = 1e-15,
    ...)
{

    ## extensive function rewriting, Giovanni Millo 29/09/2010
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

    ## change this to 'bdsmatrix'
    #require(kinship)

    ## set names for final parms vectors
    nam.beta <- dimnames(X)[[2]]
    nam.errcomp <- c("phi", "psi")

    ## initialize values for optimizer
    myparms0 <- coef0
    ## set bounds for optimizer
    lower.bounds <- c(1e-08, -0.999)
    upper.bounds <- c(1e+09, 0.999)

    ## here observations are reordered like in standard panels, to exploit
    ## the fact that in this case the vcov matrix is block-diagonal

    ## lag y parm kept for compatibility
    wy <- NULL

    ## the sigma matrix is inverted during the GLS step and not before as
    ## in the other cases, to take advantage of specialized methods in the
    ## 'kinship' (migrate to --> 'bdsmatrix'!) package

    ## GLS step function for bdsmatrices
    GLSstepBDS <- function(X, y, sigma) {
        b.hat <- solve(crossprod(X, solve(sigma, X)), crossprod(X,
            solve(sigma, y)))
        ehat <- y - X %*% b.hat
        sigma2ehat <- crossprod(ehat, solve(sigma, ehat))/(n * t)
        return(list(betahat=b.hat, ehat=ehat, sigma2=sigma2ehat))
    }

    ## rearranging module
    ## save this for eventually re-rearranging output
    oo.0 <- order(tind, ind)
    ## reorder as stacked time series, as in std. panels
    oo <- order(ind, tind)
    X <- X[oo, ]
    y <- y[oo]
    wy <- wy[oo]
    ind <- ind[oo]
    tind <- tind[oo]

    ## modules for likelihood
    Vmat <- function(rho, t) {
        V1 <- matrix(ncol = t, nrow = t)
        for (i in 1:t) V1[i, ] <- rho^abs(1:t - i)
        V <- (1/(1 - rho^2)) * V1
    }
    alfa2 <- function(rho) (1 + rho)/(1 - rho)
    d2 <- function(rho, t) alfa2(rho) + t - 1
    bSigma <- function(phirho, n, t) {
        phi <- phirho[1]
        rho <- phirho[2]
        Jt <- matrix(1, ncol = t, nrow = t)
        bSigma <- phi * Jt + Vmat(rho, t)
        bSigma
    }
    fullSigma <- function(phirho, n, t) {
        phi <- phirho[1]
        rho <- phirho[2]
        ## psi not used: here passing 3 parms, but works anyway
        ## because psi is last one
        sigma.i <- bSigma(phirho, n, t)
        fullSigma <- bdsmatrix(rep(t, n), rep(as.numeric(sigma.i),
            n))
        fullSigma
    }


    ## likelihood function, both steps included
    ll.c <- function(phirho, y, X, n, t, w, w2, wy) {
        ## retrieve parms
        phi <- phirho[1]
        rho <- phirho[2]
        ## calc sigma (here not inverted)
        sigma <- fullSigma(phirho, n, t)
        ## do GLS step to get e, s2e
        glsres <- GLSstepBDS(X, y, sigma)
        e <- glsres[["ehat"]]
        s2e <- glsres[["sigma2"]]
        ## calc ll
        uno <- n/2 * log(1 - rho^2)
        due <- -n/2 * log(d2(rho, t) * (1 - rho)^2 * phi + 1)
        tre <- -(n * t)/2 * log(s2e)
        cinque <- -1/(2 * s2e) * crossprod(e, solve(sigma, e))
        const <- -(n * t)/2 * log(2 * pi)
        ll.c <- const + uno + due + tre + cinque
        ## invert sign for minimization
        llc <- -ll.c
    }

    ## generic-ssr from here

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
    sigma <- fullSigma(myparms, n, t)
    beta <- GLSstepBDS(X, y, sigma)

    ## final vcov(beta)
    covB <- as.numeric(beta[[3]]) *
        solve(crossprod(X, solve(sigma, X)))
    ## final vcov(errcomp)
    covTheta <- solve(-fdHess(myparms, function(x) -ll.c(x,
        y, X, n, t, w, w2, wy))$Hessian)          # lag-specific line: wy
    covAR <- NULL
    covPRL <- covTheta

    ## final parms
    betas <- as.vector(beta[[1]])
    sigma2 <- as.numeric(beta[["sigma2"]])    
    arcoef <- NULL
    errcomp <- myparms
    names(betas) <- nam.beta
    names(errcomp) <- nam.errcomp

    dimnames(covB) <- list(nam.beta, nam.beta)
    dimnames(covPRL) <- list(names(errcomp), names(errcomp))

    ## remember to rearrange any output as x <- x[oo.0]

    ## result
    RES <- list(betas = betas, arcoef=arcoef, errcomp = errcomp,
                covB = covB, covAR=covAR, covPRL = covPRL, ll = myll,
                sigma2 = sigma2)

    return(RES)
}
