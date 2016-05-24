semsrmod <-
function (X, y, ind, tind, n, k, t., nT, w, w2, coef0 = rep(0, 3),
    hess = FALSE, trace = trace, x.tol = 1.5e-18, rel.tol = 1e-15,
    method="nlminb", ...)
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

    ## if w2!=w has been specified, then let w=w2
    w <- w2

    ## set names for final parms vectors
    nam.beta <- dimnames(X)[[2]]
    nam.errcomp <- c("psi", "rho")

    ## initialize values for optimizer
    myparms0 <- coef0
    ## set bounds for optimizer

    ## modules for likelihood
    Vmat.1 <- function(rho, t.) {
        ## V^(-1) is 'similar' to its 3x3 counterpart,
        ## irrespective of t.:
        ## see Vmat.R in /sparsealgebra
        if(t.==1) {Vmat.1 <- 1} else {
            Vmat.1 <- matrix(0, ncol = t., nrow = t.)
            ## non-extreme diag. elements
            for (i in 2:(t.-1)) Vmat.1[i,i] <- (1-rho^4)/(1-rho^2)
            ## extremes of diagonal
            Vmat.1[1,1] <- Vmat.1[t.,t.] <- 1
            ## bidiagonal elements
            for (j in 1:(t.-1)) Vmat.1[j+1,j] <- -rho
            for (k in 1:(t.-1)) Vmat.1[k,k+1] <- -rho
        }
        return(Vmat.1)
    }
    BB.1 <- function(lambda, w) {
        solve(xprodB(lambda, listw=w))
    }
    invSigma <- function(rholambda, n, t., w) {
        ## retrieve parms
        rho <- rholambda[1]
        lambda <- rholambda[2]
        ## psi not used: here passing 4 parms, but works anyway
        ## because psi is last one
        ## calc inverse
        invVmat <- Vmat.1(rho, t.)
        BB <- xprodB(lambda, w)
        invSigma <- kronecker(invVmat, BB)
        invSigma
    }
    ## likelihood function, both steps included
    ll.c <- function(rholambda, y, X, n, t., w, w2, wy) {
        ## retrieve parms
        rho <- rholambda[1]
        lambda <- rholambda[2]
        ## calc inverse sigma
        sigma.1 <- invSigma(rholambda, n, t., w2)
        ## lag y
        ## do GLS step to get e, s2e
        glsres <- GLSstep(X, y, sigma.1)
        e <- glsres[["ehat"]]
        s2e <- glsres[["sigma2"]]
        ## calc ll
        zero <- 0
        uno <- n/2 * log(1 - rho^2)
        tre <- -(n * t.)/2 * log(s2e)
        quattro <- t. * ldetB(lambda, w2)
        cinque <- -1/(2 * s2e) * t(e) %*% sigma.1 %*% e
        const <- -(n * t.)/2 * log(2 * pi)
        ll.c <- const + zero + uno + tre + quattro + cinque
        ## invert sign for minimization
        llc <- -ll.c
    }

    ## set bounds for optimizer
    lower.bounds <- c(-0.999, -0.999)
    upper.bounds <- c(0.999, 0.999)

    ## constraints as cA %*% theta + cB >= 0
    ## equivalent to: phi>=0, -1<=(rho, lambda, psi)<=1
    ## NB in maxLik() optimization cannot start at the boundary of the
    ## parameter space !
    cA <- cbind(c(1,-1,rep(0,2)),
               c(rep(0,2), 1, -1))
    cB <- rep(1,4)

    ## generic from here

    ## GLS step function
    GLSstep <- function(X, y, sigma.1) {
        b.hat <- solve(t(X) %*% sigma.1 %*% X,
                       t(X) %*% sigma.1 %*% y)
        ehat <- y - X %*% b.hat
        sigma2ehat <- (t(ehat) %*% sigma.1 %*% ehat)/(n * t.)
        return(list(betahat=b.hat, ehat=ehat, sigma2=sigma2ehat))
    }

    ## optimization

    ## adaptive scaling
    parscale <- 1/max(myparms0, 0.1)

    if(method=="nlminb") {

        optimum <- nlminb(start = myparms0, objective = ll.c,
                          gradient = NULL, hessian = NULL,
                          y = y, X = X, n = n, t. = t., w = w, w2 = w2,
                          scale = parscale,
                          control = list(x.tol = x.tol,
                                 rel.tol = rel.tol, trace = trace),
                          lower = lower.bounds, upper = upper.bounds)

        ## log likelihood at optimum (notice inverted sign)
        myll <- -optimum$objective
        ## retrieve optimal parms and H
        myparms <- optimum$par
        myHessian <- fdHess(myparms, function(x) -ll.c(x,
                            y, X, n, t., w, w2))$Hessian

    } else {

        #require(maxLik)

        ## initial values are not allowed to be zero
        maxout<-function(x,a) ifelse(x>a, x, a)
        myparms0 <- maxout(myparms0, 0.01)

        ## invert sign for MAXimization
        ll.c2 <- function(phirholambda, y, X, n, t., w, w2) {
            -ll.c(phirholambda, y, X, n, t., w, w2)
        }

        ## max likelihood
        optimum <- maxLik(logLik = ll.c2,
                          grad = NULL, hess = NULL, start=myparms0,
                          method = method,
                          parscale = parscale,
                          constraints=list(ineqA=cA, ineqB=cB),
                          y = y, X = X, n = n, t. = t., w = w, w2 = w2)

        ## log likelihood at optimum (notice inverted sign)
        myll <- optimum$maximum  # this one MAXimizes
        ## retrieve optimal parms and H
        myparms <- optimum$estimate
        myHessian <- optimum$hessian
    }

    ## one last GLS step at optimal vcov parms
    sigma.1 <- invSigma(myparms, n, t., w)
    beta <- GLSstep(X, y, sigma.1)

    ## final vcov(beta)
    covB <- as.numeric(beta[[3]]) *
        solve(t(X) %*% sigma.1 %*% X)

    ## final vcov(errcomp)
    nvcovpms <- length(nam.errcomp) - 1
    ## error handler here for singular Hessian cases
    covTheta <- try(solve(-myHessian), silent=TRUE)
    if(class(covTheta) == "try-error") {
        covTheta <- matrix(NA, ncol=nvcovpms+1,
                           nrow=nvcovpms+1)
        warning("Hessian matrix is not invertible")
    }
    covAR <- NULL
    covPRL <- covTheta

    ## final parms
    betas <- as.vector(beta[[1]])
    sigma2 <- as.numeric(beta[["sigma2"]])
    arcoef <- NULL
    errcomp <- myparms[which(nam.errcomp!="lambda")]
    names(betas) <- nam.beta
    names(errcomp) <- nam.errcomp[which(nam.errcomp!="lambda")]

    dimnames(covB) <- list(nam.beta, nam.beta)
    dimnames(covPRL) <- list(names(errcomp), names(errcomp))

    ## result
    RES <- list(betas = betas, arcoef=arcoef, errcomp = errcomp,
                covB = covB, covAR=covAR, covPRL = covPRL, ll = myll,
                sigma2 = sigma2)

    return(RES)
}


