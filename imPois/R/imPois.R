#' @title Imprecise Inferential Framework for Poisson Sampling Model
#' 
#' @description  Imprecise probabilities introduced by Peter Walley in 1991 is the basis of this imprecise inferential framework.  The package provides a collection of tools for conducting analysis of epistemic uncerstainty with Poisson and zero-truncated Poisson sampling parameter estimation. 
#'
#' @name imPois
#' @docType package
NULL

#' @rdname kernel_of_imprecise_prob_measure
#' @title Kernel of Imprecise Probability Measure Formulated By Bickis and Lee 
#' @description Imprecise probability density function proposed by Bickis and Lee (2014) is defined.  See \sQuote{Details}.
#' @param t random variable
#' @param xi2 parameter associated with precision 
#' @param xi1 parameter associated with linear combination
#' @param xi0 parameter associated with effective sample size
#' @param log logical; if TRUE (default), a returned value is given in logarithm scale. 
#'
#' @details  The formal definition of Bickis and Lee's conjugate formulation is
#' \deqn{e^(-\xi_2\theta^2 + \xi_1\theta - \xi_0\exp(\theta))}
#' \eqn{\theta} is ranged from \code{-Inf} to \code{Inf}. 
#'
#' @references
#' Lee, C.H. (2014) Imprecise Prior for Imprecise Inference on Poisson Sampling Model, PhD Thesis, Biostatistics Program, University of Saskatchewan
#' 
# @return
#' 
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}
#' @export
kcpm <- function(t, xi2, xi1, xi0, log=FALSE){
	
	stopifnot(length(xi2)==1, length(xi1)==1, length(xi0)==1)
	
	val <- exp(-xi2*t^2 + xi1*t - exp(log(xi0)+t))
	if(log) val <- log(val) 
	
	return(val)
}


#' @rdname normalize_imprecise_prob_measure
#' @title Comupting Normalizing Constant of Bickis and Lee's Probability Distribution
#' @description Given parameters, a normalizing constant for Bickis and Lee's probability distribution is computed. 
#' @param xi2 parameter associated with precision 
#' @param xi1 parameter associated with linear combination
#' @param xi0 effective sample size
#' @param log a logical value; if TRUE, probabilities are given as \eqn{log(p)}.
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}
#' @export 
cgf <- function(xi2, xi1, xi0, log=TRUE){
	
	stopifnot(length(xi2)==1, length(xi1)==1, length(xi0)==1)
	
	op <- obj <- list()
	cutoff <- c(Inf, 1e9, 1e8, 1e7, 5e6, 1e6, 5e5, 1e5, 5e4, 1e4, 5e3, 1e3, 5e2, 3e2, 1e2, 5e1, 3e1, 2e1)
	
	for(i in 1:length(cutoff)){
		limit <- cutoff[i]
		op[[i]] <- tryCatch(stats::integrate(kcpm, lower=-limit, upper=limit, xi2=xi2, xi1=xi1, xi0=xi0), 
			error=function(e){
				errmsg <- e$message
				if (grepl("non-finite function", errmsg)) e$value <- e$abs.error <- Inf
				else print(errmsg)
				return(e)
			})
			
		if(is.finite(op[[1]]$value)){
			value <- op[[1]]$value
			value <- if(log) log(value) else value
			return(list(attempts=op, limit=limit, value=value, i=i))
		}
			
		if (i > 2) {
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
	
	vals <- do.call(c,lapply(op, "[[", "value"))
	value <- if(log) log(vals[i]) else vals[i]
	names(vals) <- 1:i
	
	robj <- list(vals=vals, attempts=op, limit=limit, value=value, i=(i-2))
	return(robj)
}


#' @rdname imprecise_distribution
#' @title Imprecise Probability Distribution
#' @description Density and distribution function for the Bickis and Lee's distribution with three parameters \eqn{xi_2}, \eqn{xi_1}, and \eqn{xi_0}. 
#' @param x quantiles
#' @param pars a numeric vector of parameters
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca}
#' @export 
dcpm <- function(x, pars){
	
	stopifnot(length(pars) == 3)
	
	p2 <- pars[1]
	p1 <- pars[2]
	p0 <- pars[3]
	fx <- kcpm(t=x, xi2=p2, xi1=p1, xi0=p0)
	nconst <- cgf(xi2=p2, xi1=p1, xi0=p0)$value
	if( all(is.finite(fx), is.finite(nconst))) p <- fx/exp(nconst) else p <- 0
	return(p)
}

#' @rdname imprecise_distribution
#' @param q quantiles
#' @export
pcpm <- function(q, pars){
	
	stopifnot(length(pars)==3)
	
	rv <- stats::integrate(dcpm, lower=-Inf, upper=q, pars=pars)$value
	
	return(rv)
}
pcpm <- Vectorize(pcpm, "q")


#' @rdname expected_value_of_imprecise_distribution
#' @title Expected Value of Canonical Variable
#' @description Expected value of canonical variable is computed using \code{integrate} function.
#' @param y a vector of observations
#' @param pars a numeric vector of parameters 
#' @author Chel Hee Lee \email{chl948@@mail.usask.ca} 
#' @export
evfn <- function(y, pars){
	
	stopifnot(length(pars) == 3, is.vector(y))
	
	xi2 <- pars[1]
	xi1 <- pars[2]
	xi0 <- pars[3]
	
	pars1 <- c(xi2, xi1+sum(y), xi0+length(y))
	names(pars1) <- names(pars) <- c("xi2", "xi1", "xi0")
	
	op <- cgf(xi2=pars1[1], xi1=pars1[2], xi0=pars1[3])
	limit <- op$limit
	ft <- function(t, ...) t*dcpm(t, pars=pars1)
	value <- stats::integrate(ft, lower=-limit, upper=limit)$value
	
	robj <- list(y=y, pars=pars, pars1=pars1, attempts=op, value=value)
	return(robj)
}
