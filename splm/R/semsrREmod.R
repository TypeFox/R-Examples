semsrREmod <-
function (X, y, ind, tind, n, k, t., nT, w, w2, coef0 = rep(0, 4),
    hess = FALSE, trace = trace, x.tol = 1.5e-18, rel.tol = 1e-15,
    method="nlminb",
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

    ## from version 2: both nlminb (fastest, some negative covariances) and
    ## optimizers from maxLik ("BFGS", "SANN", "NM") are supported.

    ## from version 3: analytical inverse for Vmat (see notes; eliminates
    ## singular matrix pbs., from ca. 45'' to ca. 30'' on 281x3 'datiNY' example)
    ## and exploit the fact that solve(crossprod(B))=tcrossprod(solve(B))
    ## (this last not giving any benefit; but check sparse methods on B)

    ## this version 4 (5/3/2013): sparse matrix methods if w, w2 is a 'listw'
    ## needs ldetB(), solveB(), xprodB()
    ##
    ## almost no gain on medium-sized listwNY example, T=3
    
    ## if w2!=w has been specified, then let w=w2
    w <- w2 # uses w2 everywhere, but just to be sure...
    
    ## set names for final parms vectors
    nam.beta <- dimnames(X)[[2]]
    nam.errcomp <- c("phi", "psi", "rho")

    ## initialize values for optimizer
    myparms0 <- coef0

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
    alfa2 <- function(rho) (1 + rho)/(1 - rho)
    d2 <- function(rho, t.) alfa2(rho) + t. - 1
    Jt <- matrix(1, ncol = t., nrow = t.)
    In <- diag(1, n)
    det2 <- function(phi, rho, lambda, t., w) det(d2(rho, t.) * (1 -
        rho)^2 * phi * In + BB.1(lambda, w))
    Z0 <- function(phi, rho, lambda, t., w) solve(d2(rho, t.) * (1 -
        rho)^2 * phi * In + BB.1(lambda, w))
    invSigma <- function(phirholambda, n, t., w) {
        ## retrieve parms
        phi <- phirholambda[1]
        rho <- phirholambda[2]
        lambda <- phirholambda[3]
        ## psi not used: here passing 4 parms, but works anyway
        ## because psi is last one
        ## calc inverse
        invVmat <- Vmat.1(rho, t.)    #
        BB <- xprodB(lambda, w)
        invSi1 <- kronecker(invVmat, BB)
        invSi2 <- 1/(d2(rho, t.) * (1 - rho)^2)
        invSi3 <- kronecker(invVmat %*% Jt %*% invVmat,  #
            Z0(phi, rho, lambda, t., w) - BB)
        invSigma <- invSi1 + invSi2 * invSi3
        invSigma
    }
    ## likelihood function, both steps included
    ll.c <- function(phirholambda, y, X, n, t., w, w2) {
        ## retrieve parms
        phi <- phirholambda[1]
        rho <- phirholambda[2]
        lambda <- phirholambda[3]
        ## calc inverse sigma
        sigma.1 <- invSigma(phirholambda, n, t., w2)
        ## do GLS step to get e, s2e
        glsres <- GLSstep(X, y, sigma.1)
        e <- glsres[["ehat"]]
        s2e <- glsres[["sigma2"]]
        ## calc ll
        zero <- 0
        uno <- n/2 * log(1 - rho^2)
        due <- -1/2 * log(det2(phi, rho, lambda, t., w2))
        tre <- -(n * t.)/2 * log(s2e)
        quattro <- (t. - 1) * ldetB(lambda, w2)
        cinque <- -1/(2 * s2e) * crossprod(e, sigma.1) %*% e
        const <- -(n * t.)/2 * log(2 * pi)
        ll.c <- const + zero + uno + due + tre + quattro + cinque
        ## invert sign for minimization
        llc <- -ll.c
    }

    ## set bounds for optimizer
    lower.bounds <- c(1e-08, -0.999, -0.999)
    upper.bounds <- c(1e+09, 0.999, 0.999)

    ## constraints as cA %*% theta + cB >= 0
    ## equivalent to: phi>=0, -1<=(rho, lambda, psi)<=1
    ## NB in maxLik() optimization cannot start at the boundary of the
    ## parameter space !
    cA <- cbind(c(1, rep(0,4)),
               c(0,1,-1,rep(0,2)),
               c(rep(0,3), 1, -1))
    cB <- c(0, rep(1,4))

    ## generic from here

    ## GLS step function
    GLSstep <- function(X, y, sigma.1) {
        b.hat <- solve(crossprod(X, sigma.1) %*% X,
                       crossprod(X, sigma.1) %*% y)
        ehat <- y - X %*% b.hat
        sigma2ehat <- (crossprod(ehat, sigma.1) %*% ehat)/(n * t.)
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
    sigma.1 <- invSigma(myparms, n, t., w2)
    beta <- GLSstep(X, y, sigma.1)

    ## final vcov(beta)
    covB <- as.numeric(beta[[3]]) *
        solve(crossprod(X, sigma.1) %*% X)

    ## final vcov(errcomp)
    nvcovpms <- length(nam.errcomp)
    ## error handler here for singular Hessian cases
    covTheta <- try(solve(-myHessian), silent=TRUE)
    if(class(covTheta) == "try-error") {
        covTheta <- matrix(NA, ncol=nvcovpms,
                           nrow=nvcovpms)
        warning("Hessian matrix is not invertible")
    }
    covAR <- NULL
    covPRL <- covTheta[1:nvcovpms, 1:nvcovpms, drop=FALSE]

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
