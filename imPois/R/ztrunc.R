#' @rdname kernel_of_imprecise_prob_measure
#' @param m random variable
#' @export
kcpm_m <- function(m, xi2, xi1, xi0, log=FALSE){
	val <- exp(-xi2*log(m)^2 + (xi1-1)*log(m) - exp(log(xi0)+log(m)))
	if(log) val <- log(val) 
	return(val)
}

#' @rdname kernel_of_imprecise_prob_measure
#' @param ny number of observations
#' @export 
kcpm.ztrunc <- function(m, xi2, xi1, xi0, ny, log=FALSE){
	lv <- kcpm_m(m=m, xi2=xi2, xi1=xi1, xi0=xi0, log=TRUE)
	ln.1mp0 <- stats::ppois(q=0, lambda=m, lower.tail=FALSE, log.p=TRUE)
	lv <- lv - ny*ln.1mp0
	val <- exp(lv)
	if(log) val <- log(val)
	return(val)
}

#' @rdname normalize_imprecise_prob_measure
#' @param ny number of observations
#' @export 
cgf.ztrunc <- function(xi2, xi1, xi0, ny){ 
	ev <- stats::integrate(kcpm.ztrunc, lower=0, upper=Inf, xi2=xi2, xi1=xi1, xi0=xi0, ny=ny)
	ev <- log(ev$value)
	robj <- list(value=ev)
	return(robj)
} 

#' @rdname imprecise_distribution
#' @param m random variable 
#' @param ny number of observations
#' @export 
dcpm.ztrunc <- function(m, pars, ny){
	fx <- kcpm.ztrunc(m=m, xi2=pars[1], xi1=pars[2], xi0=pars[3], ny=ny)
	nconst <- cgf.ztrunc(xi2=pars[1], xi1=pars[2], xi0=pars[3], ny=ny)$value
	if (all(is.finite(fx), is.finite(nconst))) p <- fx/exp(nconst) else p <- 0
	return(p)
}

#' @rdname imprecise_distribution
#' @export
pcpm.ztrunc <- function(q, pars, ny){
	
	f.ztrunc.t <- function(t, xi2, xi1, xi0, ny, log=FALSE){
		lv <- kcpm(t=t, xi2=xi2, xi1=xi1, xi0=xi0, log=TRUE)
		lv <- lv - ny*stats::ppois(q=0, lambda=exp(t), lower.tail=FALSE, log.p=TRUE)
		val <- exp(lv)
		return(val)
	}
	dft <- function(t, pars, ny){
		ft <- f.ztrunc.t(t=t, xi2=pars[1], xi1=pars[2], xi0=pars[3], ny=ny, log=FALSE)
		nconst <- cgf.ztrunc(xi2=pars[1],xi1=pars[2],xi0=pars[3],ny=ny)$value
		if (all(is.finite(ft), is.finite(nconst))) p <- ft/exp(nconst) else p <- 0
		return(p)
	}

	err <- FALSE 
	i <- 1
	cutoff <- c(Inf, 1e4, 5e3, 3e3, 1e3, 9e2, 8e2, 7e2, 6e2, 5e2, 4e2, 3e2, 2e2, 1e2)
	op <- list()
	for (i in 1:length(cutoff)) {
		limit <- cutoff[i]
		op[[i]] <- tryCatch(stats::integrate(dft, lower=-limit, upper=q, pars=pars, ny=ny), error=function(e) {
			e$value <- Inf
			return(e)
			})
		if (i > 2){
			v_2 <- op[[i-2]]$value
			v_1 <- op[[i-1]]$value
			v0 <- op[[i]]$value
			err01 <- abs(v0-v_1)
			err02 <- abs(v0-v_2)
			err12 <- abs(v_1-v_2)
			if(is.na(err01)) next
			if ( (err01<1e-8) & (err02<1e-8)) break
		}
	}
	rv <- op[[i]]$value
	return(rv)
}
pcpm.ztrunc <- Vectorize(pcpm.ztrunc, "q")

#' @rdname expected_value_of_imprecise_distribution
#' @export 
evfn.ztrunc <- function(y, pars){
	ny <- length(y)
	pars1 <- c(pars[1], pars[2]+sum(y), pars[3]+ny)
	names(pars1) <- names(pars) <- c("xi2", "xi1", "xi0")
	
	ft <- function(m, ...) log(m)*dcpm.ztrunc(m, pars=pars1, ny=ny)
	value <- stats::integrate(ft, lower=0, upper=Inf)$value
	
	robj <- list(y=y, pars=pars, pars1=pars1, value=value)
	return(robj)
}

# E(X) = lambda/(1-exp(-lambda)) = (lambda*exp(lambda))/(exp(lambda)-1)
# V(X) = E(X)[1- lambda/(exp(lambda)-1)]


#' @rdname kernel_of_imprecise_prob_measure
#' @export
kcpm.ztrunc_t <- function(t, xi2, xi1, xi0, ny, log=TRUE){
	lv <- kcpm(t=t, xi2=xi2, xi1=xi1, xi0=xi0, log=TRUE)
	lv <- lv - ny*stats::ppois(q=0, lambda=exp(t), lower.tail=FALSE, log.p=TRUE)
	return(lv)
}

#' @rdname expected_value_of_imprecise_distribution
#' @param len.chain a length of MH chain
#' @param len.burnin a length of burn-in period 
#' @export
mh.ztrunc <- function(y, pars, len.chain=1e4, len.burnin=1e3){
	xi2 <- pars[1]
	xi1 <- pars[2] + sum(y)
	xi0 <- pars[3] + length(y)
	ny <- length(y)
	len.tot <- len.chain + len.burnin
	theta <- numeric(len.tot)
	accept <- 1
	for( t in 2:len.tot){
		curr <- theta[t-1]
		cand <- stats::rnorm(1, mean=curr, sd=sqrt(0.5))
		lp.cand <- kcpm.ztrunc_t(t=cand, xi2=xi2, xi1=xi1, xi0=xi0, ny=ny, log=TRUE)
		lp.curr <- kcpm.ztrunc_t(t=curr, xi2=xi2, xi1=xi1, xi0=xi0, ny=ny, log=TRUE)
		ratio <- exp(lp.cand - lp.curr)
		if ( stats::runif(1) < ratio ) {
			theta[t] <- cand
			accept <- c(accept,1)
		} else {
			theta[t] <- curr
			accept <- c(accept,0)
		}
	}
	return(list(chain=theta, value=mean(theta[-(1:len.burnin)]), accept=accept))
}

#' @rdname expected_value_of_imprecise_distribution
#' @param const adjustment constant
#' @export 
lapprox.ztrunc <- function(y, pars, const){
	xi2 <- pars[1]
	xi1 <- pars[2] + sum(y)
	xi0 <- pars[3] + length(y)
	ny <- length(y)
	
	fn <- function(t, numer=FALSE, xi2, xi1=xi1, xi0=xi0, ny=ny,log=TRUE, ...){
		lv <- kcpm(t=t, xi2=xi2, xi1=xi1, xi0=xi0, log=TRUE)
		lv <- lv - ny*stats::ppois(q=0, lambda=exp(t), lower.tail=FALSE, log.p=TRUE)	
		if(numer) lv <- lv + log(t + const)
		return(lv)
	}
	gr <- function(t, numer=FALSE, xi2, xi1, xi0, ny, log=TRUE, ...){
		if (!numer) lv <- -2*xi2*t + xi1 - exp(log(xi0)+t) - ny*exp(t)*stats::ppois(q=0, lambda=exp(t), lower.tail=TRUE, log.p=FALSE)/stats::ppois(q=0, lambda=exp(t), lower.tail=FALSE, log.p=FALSE)
		if (numer) lv <- 1/(t + const) -2*xi2*t + xi1 - exp(log(xi0)+t) - ny*exp(t)*stats::ppois(q=0, lambda=exp(t), lower.tail=TRUE, log.p=FALSE)/stats::ppois(q=0, lambda=exp(t), lower.tail=FALSE, log.p=FALSE)
		return(lv)
	}
	fit0 <- stats::optim(par=mean(y), fn=fn, gr=gr, xi2=xi2, xi1=xi1, xi0=xi0, ny=ny, numer=FALSE, log=TRUE, method="BFGS", control=list(fnscale=-1), hessian=TRUE)
	fit1 <- stats::optim(par=mean(y), fn=fn, gr=gr, xi2=xi2, xi1=xi1, xi0=xi0, ny=ny, numer=TRUE, log=TRUE, method="BFGS", control=list(fnscale=-1), hessian=TRUE)
	sigma0 <- -1/fit0$hessian
	sigma1 <- -1/fit1$hessian
	sratio <- sigma1/sigma0
	fdiff <- exp(fit1$value - fit0$value)
	theta <- sqrt(sratio)*fdiff - const
	return(list(value=theta, fit0=fit0, fit1=fit1))
}

# E(X) = lambda/(1-exp(-lambda)) = (lambda*exp(lambda))/(exp(lambda)-1)
# V(X) = E(X)[1- lambda/(exp(lambda)-1)]
