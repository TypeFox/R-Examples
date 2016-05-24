sarsrmod <-
function (X, y, ind, tind, n, k, t, nT, w, w2, coef0 = rep(0, 2),
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

    # mark
    #print("uso versione 1") # fixed vcov.arcoef, was missing

    ## set names for final parms vectors
    nam.beta <- dimnames(X)[[2]]
    nam.errcomp <- c("psi","lambda")

    ## initialize values for optimizer
    myparms0 <- coef0
    ## set bounds for optimizer
    lower.bounds <- c(-0.999, -0.999)      # lag-specific line (2nd parm)
    upper.bounds <- c(0.999, 0.999)        # lag-specific line (idem)

    ## here first y is lagged using the data as sent from spreml()
    ## then observations are reordered like in standard panels, to exploit
    ## the fact that in this case the vcov matrix is block-diagonal

    ## calc. Wy (spatial lag of y)
    ## (flexible fun accepting either listws or matrices for w)
    Wy <- function(y, w, tind) {                  # lag-specific line
        wyt <- function(y, w) {                   # lag-specific line
            if("listw" %in% class(w)) {           # lag-specific line
                wyt <- lag.listw(w, y)            # lag-specific line
            } else {                              # lag-specific line
                wyt <- w %*% y                    # lag-specific line
            }                                     # lag-specific line
            return(wyt)                           # lag-specific line
        }                                         # lag-specific line
        wy<-list()                                # lag-specific line
        for (j in 1:length(unique(tind))) {       # lag-specific line
             yT<-y[tind==unique(tind)[j]]         # lag-specific line
             wy[[j]] <- wyt(yT, w)                # lag-specific line
             }                                    # lag-specific line
        return(unlist(wy))                        # lag-specific line
    }                                             # lag-specific line

    ## lag y once for all
    wy <- Wy(y, w, tind)                          # lag-specific line

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
        ## this is kept for standardization, as here bSigma=Vmat
        rho <- phirho[1]
        bSigma <- Vmat(rho, t) # simply
        bSigma
    }
    fullSigma <- function(phirho, n, t) {
        rho <- phirho[1]
        ## psi not used: here passing 2 parms, but works anyway
        ## because psi is last one
        sigma.i <- bSigma(phirho, n, t)
        fullSigma <- bdsmatrix(rep(t, n), rep(as.numeric(sigma.i),
            n))
        fullSigma
    }


    ## likelihood function, both steps included
    ll.c <- function(phirho, y, X, n, t, w, w2, wy) {
        ## retrieve parms
        rho <- phirho[1]
        psi <- phirho[2]                          # lag-specific line
        ## calc sigma (here not inverted)
        sigma <- fullSigma(phirho, n, t)
        ## lag y
        Ay <- y - psi * wy                        # lag-specific line
        ## do GLS step to get e, s2e
        glsres <- GLSstepBDS(X, Ay, sigma)        # lag-specific line (Ay for y)
        e <- glsres[["ehat"]]
        s2e <- glsres[["sigma2"]]
        ## calc ll
        zero <- t*ldetB(psi, w)              # lag-specific line (else zero <- 0)
        uno <- n/2 * log(1 - rho^2)
        tre <- -(n * t)/2 * log(s2e)
        cinque <- -1/(2 * s2e) * crossprod(e, solve(sigma, e))
        const <- -(n * t)/2 * log(2 * pi)
        ll.c <- const + zero + uno + tre + cinque
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
    Ay <- y - myparms[length(myparms)] * wy       # lag-specific line
    beta <- GLSstepBDS(X, Ay, sigma)

    ## final vcov(beta)
    covB <- as.numeric(beta[[3]]) *
        solve(crossprod(X, solve(sigma, X)))
    ## final vcov(errcomp)
    covTheta <- solve(-fdHess(myparms, function(x) -ll.c(x,
        y, X, n, t, w, w2, wy))$Hessian)          # lag-specific line: wy
    nvcovpms <- length(nam.errcomp) - 1
    covAR <- covTheta[nvcovpms+1, nvcovpms+1, drop=FALSE]
    covPRL <- covTheta[1:nvcovpms, 1:nvcovpms, drop=FALSE]

    ## final parms
    betas <- as.vector(beta[[1]])
    sigma2 <- as.numeric(beta[["sigma2"]])
    arcoef <- myparms[which(nam.errcomp=="lambda")]  # lag-specific line
    errcomp <- myparms[which(nam.errcomp!="lambda")]
    names(betas) <- nam.beta
    names(arcoef) <- "lambda"                        # lag-specific line
    names(errcomp) <- nam.errcomp[which(nam.errcomp!="lambda")]

    dimnames(covB) <- list(nam.beta, nam.beta)
    dimnames(covAR) <- list(names(arcoef), names(arcoef))
    dimnames(covPRL) <- list(names(errcomp), names(errcomp))

    ## remember to rearrange any output as x <- x[oo.0]

    ## result
    RES <- list(betas = betas, arcoef=arcoef, errcomp = errcomp,
                covB = covB, covAR=covAR, covPRL = covPRL, ll = myll,
                sigma2 = sigma2)

    return(RES)
}
