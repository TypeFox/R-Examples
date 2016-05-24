sarmod <-
function (X, y, ind, tind, n, k, t, nT, w, w2, coef0 = 0,
    hess = FALSE, trace = trace, x.tol = 1.5e-18, rel.tol = 1e-15,
    ...)
{

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

    ## optimized with OLS step instead of GLS

    ## set names for final parms vectors
    nam.beta <- dimnames(X)[[2]]
    nam.errcomp <- c("lambda")

    ## initialize values for optimizer
    myparms0 <- coef0
    ## set bounds for optimizer
    lower.bounds <- c(-0.999)    # lag-specific line (3rd parm)
    upper.bounds <- c(0.999)     # lag-specific line (idem)

    ## modules for likelihood
    B <- function(lambda, w) diag(1, ncol(w)) - lambda * w
    detB <- function(lambda, w) det(B(lambda, w))

    ## likelihood function, both steps included
    ll.c <- function(phipsi, y, X, n, t, w, w2, wy) {
        ## retrieve parms
        psi <- phipsi                             # lag-specific line
        ## calc inverse sigma not needed (is I_n)
        ## lag y
        Ay <- y - psi * wy                        # lag-specific line
        ## do OLS step to get e, s2e
        glsres <- OLSstep(X, Ay)                  # lag-specific line (Ay for y)
        e <- glsres[["ehat"]]
        s2e <- glsres[["sigma2"]]
        ## calc ll
        zero <- t*ldetB(psi, w2)              # lag-specific line (else zero <- 0)
        due <- 0
        tre <- -n * t/2 * log(s2e)
        quattro <- -1/(2 * s2e) * crossprod(e)
        const <- -(n * t)/2 * log(2 * pi)
        ll.c <- const + zero + due + tre + quattro
        ## invert sign for minimization
        llc <- -ll.c
    }

    ## generic from here

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

    ## OLS step function (notice: OLS here)
    OLSstep <- function(X, y) {
        b.hat <- solve(crossprod(X), crossprod(X, y))
        ehat <- y - X %*% b.hat
        sigma2ehat <- crossprod(ehat)/(n * t)
        return(list(betahat=b.hat, ehat=ehat, sigma2=sigma2ehat))
    }

    ## lag y once for all
    wy <- Wy(y, w, tind)                          # lag-specific line

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
    Ay <- y - myparms[length(myparms)] * wy       # lag-specific line
    beta <- OLSstep(X, Ay)

    ## final vcov(beta)
    covB <- as.numeric(beta[[3]]) *
        solve(crossprod(X))

    ## final vcov(errcomp)
    covTheta <- solve(-fdHess(myparms, function(x) -ll.c(x,
        y, X, n, t, w, w2, wy))$Hessian)          # lag-specific line: wy
    nvcovpms <- length(nam.errcomp) - 1
    covAR <- covTheta[nvcovpms+1, nvcovpms+1, drop=FALSE]
    covPRL <- NULL                                # ols.errors-specific

    ## final parms
    betas <- as.vector(beta[[1]])
    sigma2 <- as.numeric(beta[["sigma2"]])
    arcoef <- myparms[which(nam.errcomp=="lambda")]  # lag-specific line
    errcomp <- NULL                               # ols.errors-specific
    names(betas) <- nam.beta
    names(arcoef) <- "lambda"                        # lag-specific line
    #names(errcomp) <- nam.errcomp[which(nam.errcomp!="lambda")]

    dimnames(covB) <- list(nam.beta, nam.beta)
    dimnames(covAR) <- list(names(arcoef), names(arcoef))
    #dimnames(covPRL) <- list(names(errcomp), names(errcomp))

    ## result
    RES <- list(betas = betas, arcoef=arcoef, errcomp = errcomp,
                covB = covB, covAR=covAR, covPRL = covPRL, ll = myll,
                sigma2 = sigma2)

    return(RES)
}
