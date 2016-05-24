## This file is part of the 'agop' library.
##
## Copyright 2013 Marek Gagolewski, Anna Cena
##
## Parts of the code are taken from the 'CITAN' R package by Marek Gagolewski
##
## 'agop' is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## 'agop' is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with 'agop'. If not, see <http://www.gnu.org/licenses/>.


#' @title Parameter Estimation in the Pareto-II Distribution (MMSE)
#' 
#' @description
#' Finds the MMS estimator of the type II Pareto distribution parameters
#' using the Bayesian method (and the R code) developed by
#' Zhang and Stevens (2009).
#'
#' 
#' @param x a non-negative numeric vector
#' @return
#' a numeric vector  with the following named components:
#' \itemize{
#' \item \code{k} - estimated parameter of shape,
#' \item \code{s} - estimated parameter of scale.
#' }
#' @export
#' @family Pareto2
#' @references
#' Zhang J., Stevens M.A., A New and Efficient Estimation Method 
#' for the Generalized Pareto Distribution, Technometrics 51(3), 2009, 316-325.\cr
pareto2_estimate_mmse <- function(x)
{
   stopifnot(is.numeric(x), length(x) >= 2, is.finite(x), x >= 0)
	n <- length(x)
	x <- sort(x)

	lx <- function(b, x)
	{
		k <- -mean(log(1-b*x))
		log(b/k)+k-1
	}

	m <- 20+floor(sqrt(n))

	b <- w <- L <- 1/x[n]+(1-sqrt(m/((1:m)-0.5)))/3/x[floor(n/4+0.5)]

	for (i in 1:m) L[i] <- n*lx(b[i],x)

	for (i in 1:m) w[i]<- 1/sum(exp(L-L[i]))

	b <- sum(b*w)

	k <- 1/mean(log(1-b*x))
	s <- -1/b

	if (k<=0)  warning("estimated shape parameter <= 0")
	if (s<=0)  warning("estimated scale parameter <= 0")

	c(k=k, s=s)
}




#' @title Parameter Estimation in the Pareto-II Distribution (MLE)
#' 
#' @description
#' Finds the maximum likelihood estimator of the type II Pareto distribution's
#' shape parameter \eqn{k} and, if not given explicitly,
#'  scale parameter \eqn{s}.
#'
#' @details
#' Note that if \eqn{s} is not given, then
#' the maximum of the likelihood function may not exist
#' for some input vectors. This estimator may have large mean squared error.
#' Consider using \code{\link{pareto2_estimate_mmse}}.
#' 
#' For known \eqn{s}, the estimator is unbiased.
#' 
#' @param x a non-negative numeric vector
#' @param s a-priori known scale parameter, \eqn{s>0} or
#' \code{NA} if unknown (default)
#' @param smin lower bound for the scale parameter to look for
#' @param smax upper bound for the scale parameter to look for
#' @param tol the desired accuracy (convergence tolerance)
#' @return
#' a numeric vector  with the following named components:
#' \itemize{
#' \item \code{k} - estimated parameter of shape
#' \item \code{s} - estimated (or known, see the \code{s} argument) parameter of scale
#' }
#' or \code{c(NA, NA)} if the maximum of the likelihood function 
#' could not be found.
#' @export
#' @family Pareto2
pareto2_estimate_mle <- function(x, s=NA_real_, smin=1e-4, smax=20, tol=.Machine$double.eps^0.25)
{
   stopifnot(is.numeric(x), length(x) >= 2, is.finite(x), x >= 0)
   
   if (!identical(s, NA_real_)) {
      stopifnot(is.numeric(s), length(s) == 1, s > 0)
      return(c(k=(length(x)-1)/sum(log(1+x/s)), s=s))
   }
   
   
	n <- length(x)

	flow <- 1+sum(log(1+x/smin))/n-n/sum(1/(1+x/smin))
	fupp <- 1+sum(log(1+x/smax))/n-n/sum(1/(1+x/smax))

	if (flow*fupp >= 0)
	{
		warning("Maximum of the likelihood function could not be found")
		return(c(k=NA_real_, s=NA_real_))
	}

	s <- uniroot( function(s,x) {
		dx <- 1+x/s
		1+sum(log(dx))/n-n/sum(1/dx)
	}, c(smin, smax), x, f.lower=flow, f.upper=fupp, tol=tol)$root
	k <- (n/sum(log(1+x/s)))

	c(k=k, s=s)
}

