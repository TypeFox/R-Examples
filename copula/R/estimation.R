## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


### Estimation for nested Archimedean copulas

### initial interval/value for optimization procedures #########################

##' Compute an initial interval or value for optimization/estimation routines
##' (only a heuristic; if this fails, choose your own interval or value)
##'
##' @title Compute an initial interval or value for estimation procedures
##' @param family Archimedean family
##' @param tau.range vector containing lower and upper admissible Kendall's tau
##' @param interval logical determining if an initial interval (the default) or
##'        an initial value should be returned
##' @param u matrix of realizations following a copula
##' @param method method for obtaining initial values
##' @param warn logical indicating whether a warning message is printed (the
##'        default) if the DMLE for Gumbel is < 1 or not
##' @param ... further arguments to cor() for method="tau.mean"
##' @return initial interval or value which can be used for optimization
##' @author Marius Hofert
initOpt <- function(family, tau.range=NULL, interval=TRUE, u,
                    method=c("tau.Gumbel", "tau.mean"), warn=TRUE, ...)
{
    cop <- getAcop(family)
    if(is.null(tau.range)){
        tau.range <- switch(cop@name, # limiting (attainable) taus that can be dealt with in estimation/optimization/root-finding
                            "AMH" = { c(0, 1/3-5e-5) }, # FIXME: closer to 1, emle's mle2 fails; note: typically, Std. Error still not available and thus profile() may fail => adjust by hand
                            "Clayton" = { c(1e-8, 0.95) },
                            "Frank" = { c(1e-8, 0.94) }, # FIXME: beyond that, estimation.gof() fails for ebeta()!
                            "Gumbel" = { c(0, 0.95) },
                            "Joe" = { c(0, 0.95) },
                            stop("unsupported family for initOpt"))
    }
    if(interval) return(cop@iTau(tau.range)) # u is not required
    stopifnot(length(dim(u)) == 2L)
    method <- match.arg(method)
    ## estimate Kendall's tau
    tau.hat <- switch(method,
                      "tau.Gumbel" = {
                          x <- apply(u, 1, max)
                          theta.hat.G <- log(ncol(u))/(log(length(x))-log(sum(-log(x)))) # direct formula from edmle for Gumbel
                          if(theta.hat.G < 1){
                              if(warn) warning("initOpt: DMLE for Gumbel = ",theta.hat.G," < 1; is set to 1")
                              theta.hat.G <- 1
			  }
                          copGumbel@tau(theta.hat.G)
                      },
                      "tau.mean" = {
                          tau.hat.mat <- cor(u, method="kendall", ...) # matrix of pairwise tau()
                          mean(tau.hat.mat[upper.tri(tau.hat.mat)]) # mean of estimated taus
                      },
                      stop("wrong method for initOpt"))
    ## truncate to range if required
    cop@iTau(pmax(tau.range[1], pmin(tau.range[2], tau.hat)))
}


### Blomqvist's beta ###########################################################

##' Compute the sample version of Blomqvist's beta,
##' see, e.g., Schmid and Schmidt (2007) "Nonparametric inference on multivariate
##' versions of Blomqvist's beta and related measures of tail dependence"
##'
##' @title Sample version of Blomqvist's beta
##' @param u matrix of realizations following the copula
##' @param scaling if TRUE then the factors 2^(d-1)/(2^(d-1)-1) and
##'                2^(1-d) in Blomqvist's beta are omitted
##' @return sample version of multivariate Blomqvist beta
##' @author Marius Hofert
betan <- function(u, scaling = FALSE) {
    less.u <- u <= 0.5
    prod1 <- apply( less.u, 1, all)
    prod2 <- apply(!less.u, 1, all)
    b <- mean(prod1 + prod2)
    if(scaling) b else {T <- 2^(ncol(u)-1); (T*b - 1)/(T - 1)}
}

beta.hat <- function(u, scaling = FALSE) { .Deprecated("betan") ; betan(u, scaling) }


##' Compute the population version of Blomqvist's beta for Archimedean copulas
##'
##' @title Population version of Blomqvist's beta for Archimedean copulas
##' @param cop acopula to be estimated
##' @param theta copula parameter
##' @param d dimension
##' @param scaling if TRUE then the factors 2^(d-1)/(2^(d-1)-1) and
##'                2^(1-d) in Blomqvist's beta are omitted
##' @return population version of multivariate Blomqvist beta
##' @author Marius Hofert & Martin Maechler
beta. <- function(cop, theta, d, scaling=FALSE) {
    j <- seq_len(d)
    diags <- cop@psi(j*cop@iPsi(0.5, theta), theta) # compute diagonals
    b <- 1 + diags[d] + if(d < 30) sum((-1)^j * choose(d, j) * diags)
    else sum((-1)^j * exp(lchoose(d, j) + log(diags)))
    if(scaling) b else { T <- 2^(d-1); (T*b - 1)/(T - 1)}
}

##' Method-of-moment-like estimation of nested Archimedean copulas based on a
##' multivariate version of Blomqvist's beta
##'
##' @title Method-of-moment-like parameter estimation of nested Archimedean copulas
##'        based on Blomqvist's beta
##' @param u matrix of realizations following the copula
##' @param cop outer_nacopula to be estimated
##' @param interval bivariate vector denoting the interval where optimization takes
##'        place
##' @param ... additional parameters for safeUroot
##' @return Blomqvist beta estimator; return value of safeUroot (more or less
##'	    equal to the return value of uniroot)
##' @author Marius Hofert
ebeta <- function(u, cop, interval=initOpt(cop@copula@name), ...) {
    stopifnot(is(cop, "outer_nacopula"), is.numeric(d <- ncol(u)), d >= 2,
              max(cop@comp) == d)
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are supported")
    ## Note: We do not need the constants 2^(d-1)/(2^(d-1)-1) and 2^(1-d) here,
    ##	     since we equate the population and sample versions of Blomqvist's
    ##       beta anyway.
    b.hat <- betan(u, scaling = TRUE)
    d <- ncol(u)
    safeUroot(function(theta) {beta.(cop@copula, theta, d, scaling=TRUE) - b.hat},
              interval=interval, Sig=+1, check.conv=TRUE, ...)
}


### Kendall's tau ##############################################################

##' Sample tau checker
##'
##' @title Check sample versions of Kendall's tau
##' @param x vector of sample versions of Kendall's tau to be checked for whether
##'        they are in the range of tau of the corresponding family
##' @param family Archimedean family
##' @return checked and (if check failed) modified x
##' @author Marius Hofert
tau.checker <- function(x, family, warn=TRUE){
    eps <- 1e-8 ## "fixed" currently, see below
    tau.range <- switch(family,
                        ## limiting (attainable) taus that can be dealt with by
			## cop<family>@iTau() *and* that can be used to construct
                        ## a corresponding copula object; checked via:
                        ## eps <- 1e-8
                        ## th <- copAMH@iTau(c(0,1/3-eps)); onacopulaL("AMH",list(th[1], 1:5)); onacopulaL("AMH",list(th[2], 1:5))
                        ## th <- copClayton@iTau(c(eps,1-eps)); onacopulaL("Clayton",list(th[1], 1:5)); onacopulaL("Clayton",list(th[2], 1:5))
                        ## th <- copFrank@iTau(c(eps,1-eps)); onacopulaL("Frank",list(th[1], 1:5)); onacopulaL("Frank",list(th[2], 1:5))
                        ## th <- copGumbel@iTau(c(0,1-eps)); onacopulaL("Gumbel",list(th[1], 1:5)); onacopulaL("Gumbel",list(th[2], 1:5))
                        ## th <- copJoe@iTau(c(0,1-eps)); onacopulaL("Joe",list(th[1], 1:5)); onacopulaL("Joe",list(th[2], 1:5))
                        "AMH" = { c(0, 1/3-eps) },
                        "Clayton" = { c(eps, 1-eps) }, # copClayton@iTau(c(eps,1-eps))
                        "Frank" = { c(eps, 1-eps) }, # copFrank@iTau(c(eps,1-eps))
                        "Gumbel" = { c(0, 1-eps) }, # copGumbel@iTau(c(0,1-eps))
                        "Joe" = { c(0, 1-eps) }, # copJoe@iTau(c(0,1-eps))
                        stop("unsupported family for initOpt"))
    toosmall <- which(x < tau.range[1])
    toolarge <- which(x > tau.range[2])
    if(warn && length(toosmall)+length(toolarge) > 0){
	r <- range(x)
	if(length(x) == 1){
            warning("tau.checker: found (and adjusted) an x value out of range (x = ",
                    x,")")
	}else{
            warning("tau.checker: found (and adjusted) x values out of range (min(x) = ",
                    r[1],", max(x) = ",r[2],")")
        }
    }
    x. <- x
    x.[toosmall] <- tau.range[1]
    x.[toolarge] <- tau.range[2]
    x.
}

##' Compute pairwise estimators for nested Archimedean copulas based on Kendall's tau
##'
##' @title Pairwise estimators for nested Archimedean copulas based on Kendall's tau
##' @param u matrix of realizations following the copula
##' @param cop outer_nacopula to be estimated
##' @param method tau.mean indicates that the average of the sample versions of
##'               Kendall's tau are computed first and then theta is determined;
##'               theta.mean stands for first computing all Kendall's tau
##'               estimators and then returning the mean of these estimators
##' @param warn logical indicating whether warnings are produced (for AMH and in
##'        general for pairwise sample versions of Kendall's tau < 0) [the default]
##'        or not
##' @param ... additional arguments to cor()
##' @return averaged pairwise cor() estimators
##' @author Marius Hofert
etau <- function(u, cop, method = c("tau.mean", "theta.mean"), warn=TRUE, ...){
    stopifnot(is(cop, "outer_nacopula"), is.numeric(d <- ncol(u)), d >= 2,
              max(cop@comp) == d)
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are supported")
    tau.hat.mat <- cor(u, method="kendall",...) # matrix of pairwise tau()
    tau.hat <- tau.hat.mat[upper.tri(tau.hat.mat)] # all tau hat's
    ## define tau^{-1}
    tau_inv <- if(cop@copula@name == "AMH")
	function(tau) cop@copula@iTau(tau, check=FALSE, warn=warn) else cop@copula@iTau
    ## check and apply iTau in the appropriate way
    method <- match.arg(method)
    switch(method,
           "tau.mean" = {
	       mean.tau.hat <- mean(tau.hat) # mean of pairwise tau.hat
	       mean.tau.hat. <- tau.checker(mean.tau.hat, family=cop@copula@name,
                                            warn=warn) # check the mean
               tau_inv(mean.tau.hat.) # Kendall's tau corresponding to the mean of the sample versions of Kendall's taus
           },
           "theta.mean" = {
	       tau.hat. <- tau.checker(tau.hat, family=cop@copula@name, warn=warn) # check all values
               mean(tau_inv(tau.hat.)) # mean of the pairwise Kendall's tau estimators
           },
       {stop("wrong method")})
}


### Minimum distance estimation ################################################

##' Distances for minimum distance estimation
##'
##' @title Distances for minimum distance estimation
##' @param u matrix of realizations (ideally) following U[0,1]^(d-1) or U[0,1]^d
##' @param method distance methods available:
##'        mde.chisq.CvM  = map to a chi-square distribution (Cramer-von Mises distance)
##'        mde.chisq.KS   = map to a chi-square distribution (Kolmogorov-Smirnov distance)
##'        mde.gamma.CvM  = map to an Erlang (gamma) distribution (Cramer-von Mises distance)
##'        mde.gamma.KS   = map to an Erlang (gamma) distribution (Kolmogorov-Smirnov distance)
##' @return distance
##' @author Marius Hofert
emde.dist <- function(u, method = c("mde.chisq.CvM", "mde.chisq.KS", "mde.gamma.CvM",
                         "mde.gamma.KS")) {
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    d <- ncol(u)
    n <- nrow(u)
    method <- match.arg(method) # match argument method
    switch(method,
           "mde.chisq.CvM" = { # map to a chi-square distribution
               y <- sort(rowSums(qnorm(u)^2))
               Fvals <- pchisq(y, d)
               weights <- (2*(1:n)-1)/(2*n)
               1/(12*n) + sum((weights - Fvals)^2)
           },
           "mde.chisq.KS" = { # map to a chi-square distribution
               y <- sort(rowSums(qnorm(u)^2))
               Fvals <- pchisq(y, d)
               i <- 1:n
               max(Fvals[i]-(i-1)/n, i/n-Fvals[i])
           },
           "mde.gamma.CvM" = { # map to an Erlang distribution
               y <- sort(rowSums(-log(u)))
               Fvals <- pgamma(y, shape = d)
               weights <- (2*(1:n)-1)/(2*n)
               1/(12*n) + sum((weights - Fvals)^2)
           },
           "mde.gamma.KS" = { # map to an Erlang distribution
               y <- rowSums(-log(u))
               Fvals <- pgamma(y, shape = d)
               i <- 1:n
               max(Fvals[i]-(i-1)/n, i/n-Fvals[i])
           },
           ## Note: The distances S_n^{(B)} and S_n^{(C)} turned out to be (far)
           ##       too slow.
           stop("wrong distance method"))
}

##' Minimum distance estimation for nested Archimedean copulas
##'
##' @title Minimum distance estimation for nested Archimedean copulas
##' @param u matrix of realizations following the copula
##' @param cop outer_nacopula to be estimated
##' @param method distance methods available, see emde.dist
##' @param interval bivariate vector denoting the interval where optimization takes
##'        place
##' @param include.K logical indicating whether the last component, K, is also
##'        used or not
##' @param repara logical indicating whether the distance function is
##'        reparameterized for the optimization
##' @param ... additional parameters for optimize
##' @return minimum distance estimator; return value of optimize
##' @author Marius Hofert
emde <- function(u, cop, method = c("mde.chisq.CvM", "mde.chisq.KS", "mde.gamma.CvM",
                         "mde.gamma.KS"), interval = initOpt(cop@copula@name),
                 include.K = FALSE, repara = TRUE, ...)
{
    stopifnot(is(cop, "outer_nacopula"), is.numeric(d <- ncol(u)), d >= 2,
              max(cop@comp) == d)
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are supported")
    method <- match.arg(method) # match argument method
    distance <- function(theta) { # distance to be minimized
        cop@copula@theta <- theta
        u. <- htrafo(u, cop=cop, include.K=include.K, n.MC=0) # transform data [don't use MC here; too slow]
        emde.dist(u., method)
    }
    if(repara){
        ## reparameterization function
        rfun <- function(x, inverse=FALSE){ # reparameterization
            switch(cop@copula@name,
                   "AMH"={
                       x
                   },
                   "Clayton"={
                       if(inverse) tanpi(x/2) else atan(x)*2/pi
                   },
                   "Frank"={
                       if(inverse) tanpi(x/2) else atan(x)*2/pi
                   },
                   "Gumbel"={
                       if(inverse) 1/(1-x) else 1-1/x
                   },
                   "Joe"={
                       if(inverse) 1/(1-x) else 1-1/x
                   },
                   stop("emde: Reparameterization got unsupported family"))
        }
        ## optimize
	opt <- optimize(function(alpha) distance(rfun(alpha, inverse=TRUE)),
                        interval=rfun(interval), ...)
	opt$minimum <- rfun(opt$minimum, inverse=TRUE)
	opt
    }else{
        optimize(distance, interval=interval, ...)
    }
}


### Diagonal maximum likelihood estimation #####################################

##' Density of the diagonal of a nested Archimedean copula
##'
##' @title Diagonal density of a nested Archimedean copula
##' @param u evaluation point in [0,1]
##' @param cop outer_nacopula
##' @param log if TRUE the log-density is evaluated
##' @return density of the diagonal of cop
##' @author Marius Hofert
dDiag <- function(u, cop, log=FALSE) {
    stopifnot(is(cop, "outer_nacopula"), (d <- max(cop@comp)) >= 2)
    if(length(cop@childCops)) {
        stop("currently, only Archimedean copulas are supported")
    }
    else ## (non-nested) Archimedean :
        ## FIXME: choose one or the other (if a family has no such slot)
        ##    dDiagA(u, d=d, cop = cop@copula, log=log)
        cop@copula@dDiag(u, theta=cop@copula@theta, d=d, log=log)
}

##' @title Generic density of the diagonal of d-dim. Archimedean copula
##' @param u evaluation point in [0, 1]
##' @param d dimension
##' @param cop acopula
##' @param log if TRUE the log-density is evaluated
##' @return density of the diagonal of cop
##' @author Martin Maechler
dDiagA <- function(u, d, cop, log=FALSE) {
    stopifnot(is.finite(th <- cop@theta), d >= 2)
    ## catch the '0' case directly; needed, e.g., for AMH:
    if(any(copAMH@name == c("AMH","Frank","Gumbel","Joe")) &&
       any(i0 <- u == 0)) {
	if(log) u[i0] <- -Inf
	u[!i0] <- dDiagA(u[!i0], d=d, cop=cop, log=log)
	return(u)
    }
    if(log) {
        log(d) + cop@absdPsi(d*cop@iPsi(u, th), th, log=TRUE) +
            cop@absdiPsi(u, th, log=TRUE)
    } else {
        d * cop@absdPsi(d*cop@iPsi(u, th), th) * cop@absdiPsi(u, th)
    }
}

##' Maximum likelihood estimation based on the diagonal of a nested Archimedean copula
##'
##' @title Maximum likelihood estimation based on the diagonal of a nested Archimedean copula
##' @param u matrix of realizations following a copula
##' @param cop outer_nacopula to be estimated
##' @param interval bivariate vector denoting the interval where optimization takes
##'        place
##' @param warn logical indicating whether a warning message is printed (the
##'        default) if the DMLE for Gumbel is < 1 or not
##' @param ... additional parameters for optimize
##' @return diagonal maximum likelihood estimator; return value of optimize
##' @author Marius Hofert
edmle <- function(u, cop, interval=initOpt(cop@copula@name), warn=TRUE, ...)
{
    stopifnot(is(cop, "outer_nacopula"), is.numeric(d <- ncol(u)), d >= 2,
              max(cop@comp) == d) # dimension
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are supported")
    x <- apply(u, 1, max) # data from the diagonal
    ## explicit estimator for Gumbel
    if(cop@copula@name == "Gumbel") {
	th.G <- log(d)/(log(length(x))-log(sum(-log(x))))
	if(!is.finite(th.G) || th.G < 1) {
            if(warn) warning("edmle: DMLE for Gumbel = ",th.G,"; not in [1, Inf); is set to 1")
            th.G <- 1
	}
	list(minimum = th.G, objective = 0) # return value of the same structure as for optimize
    } else {
        ## optimize
	nlogL <- function(theta) # -log-likelihood of the diagonal
            -sum(cop@copula@dDiag(x, theta=theta, d=d, log=TRUE))
	optimize(nlogL, interval=interval, ...)
    }
}


### (Simulated) maximum likelihood estimation ##################################

##' (Simulated) maximum likelihood estimation for nested Archimedean copulas
##' -- *Fast* version (based on optimize()) called from enacopula
##'
##' @title (Simulated) maximum likelihood estimation for nested Archimedean copulas
##' @param u matrix of realizations following the copula
##' @param cop outer_nacopula to be estimated
##' @param n.MC if > 0 SMLE is applied with sample size equal to n.MC; otherwise,
##'        MLE is applied
##' @param interval bivariate vector denoting the interval where optimization takes
##'        place
##' @param ... additional parameters for optimize
##' @return (simulated) maximum likelihood estimator; return value of optimize
##' @author Marius Hofert
.emle <- function(u, cop, n.MC=0, interval=initOpt(cop@copula@name), ...)
{
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
	stop("currently, only Archimedean copulas are supported")
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    ## optimize
    mLogL <- function(theta) { # -log-likelihood
        cop@copula@theta <- theta
	-sum(.dnacopula(u, cop, n.MC=n.MC, log=TRUE))
    }
    optimize(mLogL, interval=interval, ...)
}


##' (Simulated) maximum likelihood estimation for nested Archimedean copulas
##'
##' @title (Simulated) maximum likelihood estimation for nested Archimedean copulas
##' @param u matrix of realizations following the copula
##' @param cop outer_nacopula to be estimated
##' @param n.MC if > 0 SMLE is applied with sample size equal to n.MC; otherwise,
##'        MLE is applied
##' @param optimizer optimizer used (if optimizer=NULL (or NA), then mle (instead
##'        of mle2) is used with the provided method)
##' @param method optim's method to be used (when optimizer=NULL or "optim" and
##'        in these cases method is a required argument)
##' @param interval bivariate vector denoting the interval where optimization
##'        takes place
##' @param start list containing the initial value(s) (unfortunately required by mle2)
##' @param ... additional parameters for optimize
##' @return an "mle2" object with the (simulated) maximum likelihood estimator.
##' @author Martin Maechler and Marius Hofert
##' Note: this is the *slower* version which also allows for profiling
emle <- function(u, cop, n.MC=0, optimizer="optimize", method,
		 interval=initOpt(cop@copula@name),
                 ##vvv awkward to be needed, but it is - by mle2():
                 start = list(theta=initOpt(cop@copula@name, interval=FALSE, u=u)),
                 ...)
{
    stopifnot(is(cop, "outer_nacopula"), is.numeric(d <- ncol(u)), d >= 2,
              max(cop@comp) == d)
    ## nLL <- function(theta) { # -log-likelihood
    ##	   cop@copula@theta <- theta
    ##	   -sum(.dnacopula(u, cop, n.MC=n.MC, log=TRUE))
    ## }
    if(length(cop@childCops))
	stop("currently, only Archimedean copulas are supported")
    else ## For (*non*-nested) copulas only:
	nLL <- function(theta)  # -(log-likelihood)
	    -sum(cop@copula@dacopula(u, theta, n.MC=n.MC, log=TRUE))

    ## optimization
    if(!(is.null(optimizer) || is.na(optimizer))) {
        ## stopifnot(requireNamespace("bbmle"))
	if(optimizer == "optimize")
	    bbmle::mle2(minuslogl = nLL, optimizer = "optimize",
		 lower = interval[1], upper = interval[2],
		 ##vvv awkward to be needed, but it is - by mle2():
		 start=start, ...)
	else if(optimizer == "optim") {
	    message(" optimizer = \"optim\" -- using mle2(); consider optimizer=NULL instead")
	    bbmle::mle2(minuslogl = nLL, optimizer = "optim", method = method,
		 start=start, ...)
	}
	else ## "general"
	    bbmle::mle2(minuslogl = nLL, optimizer = optimizer, start=start, ...)
    }
    else
	## use optim() .. [which uses suboptimal method for 1D, but provides Hessian]
	mle(minuslogl = nLL, method = method, start=start, ...)
}


### Estimation wrapper #########################################################

##' Computes the pseudo-observations for the given data matrix
##'
##' @title Pseudo-observations
##' @param x matrix of random variates to be converted to pseudo-observations
##' @param na.last passed to rank()
##' @param ties.method passed to rank()
##' @param lower.tail if FALSE, pseudo-observations when apply the empirical
##'        marginal survival functions are returned.
##' @return pseudo-observations (matrix of the same dimensions as x)
##' @author Marius Hofert
pobs <- function(x, na.last = "keep",
		 ## formals(rank) works in pre-2015-10-15 and newer version of rank():
		 ties.method = eval(formals(rank)$ties.method),
		 lower.tail = TRUE) {
    ties.method <- match.arg(ties.method)
    U <- apply(x, 2, rank, na.last=na.last, ties.method=ties.method) / (nrow(x)+1)
    if(lower.tail) U else 1-U
}

##' Computes different parameter estimates for a nested Archimedean copula
##'
##' @title Estimation procedures for nested Archimedean copulas
##' @param u data matrix (of pseudo-observations or from the copula "directly")
##' @param cop outer_nacopula to be estimated
##' @param method estimation method; can be
##'        "mle"             MLE
##'        "smle"            SMLE
##'        "dmle"            MLE based on the diagonal
##'        "mde.chisq.CvM"   minimum distance estimation based on the chisq distribution and CvM distance
##'        "mde.chisq.KS"    minimum distance estimation based on the chisq distribution and KS distance
##'        "mde.gamma.CvM"   minimum distance estimation based on the Erlang distribution and CvM distance
##'        "mde.gamma.KS"    minimum distance estimation based on the Erlang distribution and KS distance
##'        "tau.tau.mean"    averaged pairwise Kendall's tau estimator
##'        "tau.theta.mean"  average of Kendall's tau estimators
##'        "beta"            multivariate Blomqvist's beta estimator
##' @param n.MC if > 0 it denotes the sample size for SMLE
##' @param interval initial optimization interval for "mle", "smle", and "dmle"
##' @param xargs additional arguments for the estimation procedures
##' @param ... additional parameters for optimize
##' @return estimated value/vector according to the chosen method
##' @author Marius Hofert
enacopula <- function(u, cop, method=c("mle", "smle", "dmle", "mde.chisq.CvM",
                              "mde.chisq.KS", "mde.gamma.CvM", "mde.gamma.KS",
                              "tau.tau.mean", "tau.theta.mean", "beta"),
                      n.MC = if(method=="smle") 10000 else 0,
                      interval=initOpt(cop@copula@name),
                      xargs=list(), ...)
{

    ## setup
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot(0 <= u, u <= 1, is(cop, "outer_nacopula"), (d <- ncol(u)) >= 2,
              max(cop@comp) == d, n.MC >= 0, is.list(xargs))
    if(length(cop@childCops))
	stop("currently, only Archimedean copulas are supported")
    if(n.MC > 0 && method != "smle")
	stop("n.MC > 0  is not applicable to method '%s'", method)
    method <- match.arg(method)

    ## main part
    res <- switch(method,
                  "mle" =            do.call(.emle, c(list(u, cop,
                  interval = interval, ...), xargs)),
                  "smle" =           do.call(.emle, c(list(u, cop, n.MC = n.MC,
                  interval = interval, ...), xargs)),
                  "dmle" =           do.call(edmle, c(list(u, cop,
                  interval = interval, ...), xargs)),
                  "mde.chisq.CvM" =  do.call(emde, c(list(u, cop, "mde.chisq.CvM",
                  interval = interval, ...), xargs)),
                  "mde.chisq.KS" =   do.call(emde, c(list(u, cop, "mde.chisq.KS",
                  interval = interval, ...), xargs)),
                  "mde.gamma.CvM" =  do.call(emde, c(list(u, cop, "mde.gamma.CvM",
                  interval = interval, ...), xargs)),
                  "mde.gamma.KS" =   do.call(emde, c(list(u, cop, "mde.gamma.KS",
                  interval = interval, ...), xargs)),
                  "tau.tau.mean" =   do.call(etau, c(list(u, cop, "tau.mean", ...),
                  xargs)),
                  "tau.theta.mean" = do.call(etau, c(list(u, cop, "theta.mean", ...),
                  xargs)),
                  "beta" =           do.call(ebeta, c(list(u, cop,
                  interval = interval, ...), xargs)),
                  stop("wrong estimation method for enacopula"))

    ## FIXME: deal with result, check details, give warnings

    ## return the estimate
    switch(method,
           "mle" =            res$minimum,
           "smle" =           res$minimum,
           "dmle" =           res$minimum,
           "mde.chisq.CvM" =  res$minimum,
           "mde.chisq.KS" =   res$minimum,
           "mde.gamma.CvM" =  res$minimum,
           "mde.gamma.KS" =   res$minimum,
           "tau.tau.mean" =   res,
           "tau.theta.mean" = res,
           "beta" =           res$root,
           stop("wrong estimation method"))

}
