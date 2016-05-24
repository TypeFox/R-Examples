##
##  Copyright (C) 2010 Friedrich Leisch
##  Copyright (C) 2013 Martin Maechler
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  A copy of the GNU General Public License is available at
##  http://www.r-project.org/Licenses/


estep.nm <- function(x, obj, par)
{
    ## Purpose: 1 E-step for data 'x' and  obj = <norMix> object
    ##                               *or*  par = our parametrization
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 17 Dec 2007, 18:53
    if(missing(obj))
        obj <- par2norMix(par)
    else {
        if(!missing(par))
            stop("cannot specify both 'par' and 'obj'")
        if(!is.norMix(obj))
            stop("'obj' must be a 'Normal Mixture' object!")
    }
    sd <- obj[,"sigma"]
    m <- length(sd)
    n <- length(x)
    xxl <- rep(x, each = m) # so we can "recycle" (mu,sig) in next line
    ## fx are the density components apart from the common factor  1/sqrt(2*pi)
    fx <- matrix(obj[,"w"] * exp(-.5 * ((xxl - obj[,"mu"]) / sd)^2) / sd,
                 n,m, byrow=TRUE)
    ## z_{i,j} =  \pi_j f_j(x_i) / sum_j' \pi_j' f_j'(x_i)
    fx / rowSums(fx)
}

mstep.nm <- function(x, z)
{
    ## Purpose: 1 M-step for data 'x' and matrix 'z'
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 18 Dec 2007, 08:45
    stopifnot(is.numeric(x), is.numeric(z), is.matrix(z),
              (n <- length(x)) > 1, n == nrow(z), ncol(z) >= 1)
    n.j <- colSums(z)
    mu   <- colSums(z * x) / n.j
    sig2 <- colSums(z* outer(x, mu, "-")^2) / n.j
    list(w = n.j / n, mu = mu, sigma = sqrt(sig2))
}

emstep.nm <- function(x, obj)
{
    ## Purpose:   E-step + M-step  for 1-dim normal mixture <norMix>
    ##    return: an 'updated' <norMix>
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  1 Jan 2008, 17:38

    if(!is.norMix(obj))
        stop("'obj' must be a 'Normal Mixture' object!")

    sd <- obj[,"sigma"]
    m <- length(sd)
    n <- length(x)

    ## E-step:  z_{i,j} :=  Prob( x_i in C_j | theta ) --------
    xxl <- rep(x, each = m) # so we can "recycle" (mu,sig) in next line
    ## fx are the density components apart from the common factor  1/sqrt(2*pi)
    fx <- matrix(obj[,"w"]/sd * exp((xxl - obj[,"mu"])^2 / (-2*sd*sd)), n,m, byrow=TRUE)
    ## z_{i,j} =  \pi_j f_j(x_i) / sum_j' \pi_j' f_j'(x_i)
    z <- fx / rowSums(fx)

    ## M-step for data 'x' and matrix 'z' --------------------
    n.j  <- colSums(z)
    mu   <- colSums(z * x) / n.j # and use the *updated* mu:
    sig2 <- colSums(z* matrix((xxl - mu)^2, n,m, byrow=TRUE)) / n.j

    norMix(w = n.j / n, mu = mu, sigma = sqrt(sig2), name="")
}
##--end---------------------------- [em]step.nm() ---------------------------

norMixEM <- function(x, m, name = NULL, sd.min = 1e-7 * diff(range(x))/m,
                     maxiter = 100, tol = sqrt(.Machine$double.eps), trace = 1)
{
    call. <- match.call()
    if(is.null(name))
	name <- sub("\\(x = ","(",
		    sub(", trace = [^,)]+", '', deparse(call.)))
    stopifnot(is.numeric(x), (n <- length(x)) >= 1,
              (maxiter <- as.integer(maxiter)[1]) >= 1)
    m <- as.integer(m)
    stopifnot((mm <- max(m)) >= 1)
    if(length(m) > 1) {
	init <- rep(m, length = n)
	m <- mm
    }
    else {
        q <- quantile(x, seq(0, 1, by = 1/m))
        init <- as.integer(cut(x, q, labels=FALSE, include.lowest=TRUE))
    }
    stopifnot(is.numeric(sd.min), length(sd.min) == 1, sd.min >= 0)

    xxl <- rep(x, each = m) # so we can "recycle"  (mu,sig)  in E-step
    ## z := the "posterior" probabilities (in each step: prior -> post. -> prior -> ..)
    z <- matrix(0, nrow = n, ncol = m)
    z[cbind(1:n, init)] <- 1

    ## offset to get log-likelihood  =      n * log(  <common factor>  )
    l.fac <- -n*(log(n) + log(2*pi)/2) ## = n * log( 1/ (n * sqrt(2pi)))
    llh <- -Inf ; conv <- FALSE
    if(trace >= 2)
        cat(sprintf("%4s: %18s %9s\n", "it.", "log-likelihood", "rel.chng"))
    for(it in 1:maxiter) {

        ## M-step for data 'x' and matrix 'z' --------------------
        n.j <- colSums(z) # n.j/n =: \hat{\pi_j}
        mu  <- colSums(z * x) / n.j
        xm2 <- (xxl - mu)^2 ## ==> sd based on *new* mu: typically faster, but *not* classical EM
        sd  <- pmax(sd.min, sqrt(colSums(z* matrix(xm2, n,m, byrow=TRUE)) / n.j))

        ## E-step:  z_{i,j} :=  Prob( x_i in C_j | theta ) --------
        ##                   = the "posterior" probabilities
        ## fx: the density components apart from the common factor  1/n/sqrt(2*pi)
        fx <- matrix(n.j/sd * exp(xm2 / (-2* sd^2)), n,m, byrow=TRUE)
        ## z_{i,j} =  \pi_j f_j(x_i) / sum_j' \pi_j' f_j'(x_i)
        z <- fx / (pr.x <- rowSums(fx))

        ## Compute log-Likelihood and check for convergence
        llh.old <- llh
        llh <- sum(log(pr.x))
        relEr <- abs(llh-llh.old)/(abs(llh) +1e-7)
	if(trace >= 2) cat(sprintf("%4d: %#18.14g %9.3e\n",
				   it, llh+l.fac, if(it > 1) relEr else NA))
        if(conv <- relEr <= tol)# '<=' : so can use '0'
            break
    }
    if(!conv)
	warning(gettextf("EM did not converge in %d iterations (with 'tol'=%g)",
			 maxiter, tol))
    else if(trace == 1)
	cat(gettextf("EM converged in %d iterations (with relative change=%g)",
		     it, relEr), "\n", sep="")
    if(min(sd) - sd.min < 1e-5*sd.min)
	warning(gettextf("some 'sd' ended up very close to 'sd.min'=%g",
			 sd.min))

    structure(norMix(mu = mu, sigma = sd, w = n.j/n, name = name),
              loglik = llh+l.fac, iter = it, tol = relEr, converged = conv,
              class = c("fitEM", "nMfit", "norMix"))
}


norMixMLE <- function(x, m, name = NULL, ## sd.min = 1e-7 * diff(range(x))/m,
                     maxiter = 100, tol = sqrt(.Machine$double.eps), trace = 2)
{
    call. <- match.call()
    if(is.null(name))
	name <- sub("\\(x = ","(",
		    sub(", trace = [^,)]+", '', deparse(call.)))
    stopifnot(is.numeric(x), (n <- length(x)) >= 1,
              (maxiter <- as.integer(maxiter)[1]) >= 1)
    m <- as.integer(m)
    stopifnot((mm <- max(m)) >= 1)
    if(length(m) > 1) {
	init <- rep(m, length = n)
	m <- mm
    }
    else {
        q <- quantile(x, seq(0, 1, by = 1/m))
        init <- as.integer(cut(x, q, labels=FALSE, include.lowest=TRUE))
    }
    ## stopifnot(is.numeric(sd.min), length(sd.min) == 1, sd.min >= 0)

    ## z := the "posterior" probabilities, here 0/1 from the grouping
    z <- matrix(0, nrow = n, ncol = m)
    z[cbind(1:n, init)] <- 1

    ## 1 M-Step, to get initial
    n.j <- colSums(z)
    mu   <- colSums(z * x) / n.j
    sig2 <- colSums(z* outer(x, mu, "-")^2) / n.j
    par <- .nM2par(mu = mu, sigma = sqrt(sig2), w = n.j / n)

    neglogl <- function(par) -llnorMix(par, x=x)

## FIXME ??  Use  mle() and stats4 --> many methods !!

    optr <- optim(par, neglogl, method = "BFGS", control =
                list(maxit=maxiter, reltol = tol,
                     trace=(trace > 0), REPORT= pmax(1, 10 %/% trace)))
    conv <- optr$convergence == 0
    it <- optr$counts[2:1] ## "gradient" (= high-level it.) first
    m.s.w <- .par2nM(optr $ par)

    if(!conv)
	warning(gettextf("MLE did not converge in %d iterations (with 'tol'=%g)",
			 maxiter, tol))
    else if(trace == 1)
	cat(gettextf("MLE converged in %d iterations", it[[1]]), "\n", sep="")
    ## sd <- m.s.w$sd
    ## if(min(sd) - sd.min < 1e-5*sd.min)
    ##     warning(gettextf("some 'sd' ended up very close to 'sd.min'=%g",
    ##     		 sd.min))
    structure(norMix(mu = m.s.w$mu, sigma = m.s.w$sd, w = m.s.w$w, name = name),
              loglik = -optr $ value, iter = it, tol = NA, converged = conv,
              class = c("fitMLE", "nMfit", "norMix"))
}
