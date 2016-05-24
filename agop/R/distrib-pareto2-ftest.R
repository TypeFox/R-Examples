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


#' @title Two-Sample F-test For Equality of Shape Parameters for Type II-Pareto Distributions
#' 
#' @description
#' Performs F-test for equality of shape parameters
#' of two samples from the Pareto type-II distributions with known
#' and equal scale parameters, \eqn{s>0}.
#'
#' @details
#' Given two samples \eqn{(X_1,...,X_n)} i.i.d. \eqn{P2(k_x,s)}
#' and \eqn{(Y_1,...,Y_m)} i.i.d. \eqn{P2(k_y,s)}
#' this test verifies the null hypothesis
#' \eqn{H_0: k_x=k_y}
#' against two-sided or one-sided alternatives, depending
#' on the value of \code{alternative}.
#' It bases on test statistic
#' \eqn{T(X,Y)=\frac{n\sum_{i=1}^m\log(1+Y_i/m)}{m\sum_{i=1}^n\log(1+X_i/n)}}{T=n/m*sum(log(1+Y/m))/sum(log(1+X/n))}
#' which, under \eqn{H_0}, has the Snedecor's F distribution with \eqn{(2m, 2n)}
#' degrees of freedom.
#'
#' Note that for \eqn{k_x < k_y}, then \eqn{X} dominates \eqn{Y} stochastically.
#'
#' 
#' @param x a non-negative numeric vector
#' @param y a non-negative numeric vector
#' @param s the known scale parameter, \eqn{s>0}
#' @param alternative indicates the alternative hypothesis and must be one of 
#' \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}
#' @param significance significance level, \eqn{0<}\code{significance}\eqn{<1} 
#' or \code{NULL}. See the Value section for details
#' @return
#' If \code{significance} is not \code{NULL}, then
#' the list of class \code{power.htest} with the following components is passed as a result:
#' \itemize{
#' \item \code{statistic} -	the value of the test statistic.
#' \item \code{result} -	either FALSE (accept null hypothesis) or TRUE (reject).
#' \item \code{alternative} -	a character string describing the alternative hypothesis.
#' \item \code{method} -	a character string indicating what type of test was performed.
#' \item \code{data.name} -	a character string giving the name(s) of the data.
#' }
#' 
#' 
#' Otherwise, the list of class \code{htest} with the following components is passed as a result:
#' \itemize{
#' \item \code{statistic} 	the value of the test statistic.
#' \item \code{p.value} 	the p-value of the test.
#' \item \code{alternative} 	a character string describing the alternative hypothesis.
#' \item \code{method} 	a character string indicating what type of test was performed.
#' \item \code{data.name} 	a character string giving the name(s) of the data.
#' }
#' @export
#' @family Pareto2
pareto2_test_f <- function(x, y, s, alternative = c("two.sided", "less", "greater"), significance=NULL)
{
	alternative <- match.arg(alternative)
	DNAME <- deparse(substitute(x))
	DNAME <- paste(DNAME, "and", deparse(substitute(y)))

	x <- x[!is.na(x)]
	nx <- length(x)
	if (nx < 1L || any(x<0)) stop("incorrect 'x' data")

	y <- y[!is.na(y)]
	ny <- length(y)
	if (ny < 1L || any(y<0)) stop("incorrect 'y' data")


	if (mode(s) != "numeric" || length(s) != 1 || s <= 0) stop("'s' should be > 0")




	STATISTIC <- nx/ny*sum(log(1+y/s))/sum(log(1+x/s))
	names(STATISTIC) <- "F"
	METHOD <- "Two-sample F-test for equality of shape parameters for Type II-Pareto distributions with known common scale parameter";
	nm_alternative <- switch(alternative, two.sided = "two-sided",
		less = "kx < ky",
		greater = "kx > ky")



	if (!is.null(significance))
	{
		if (length(significance) != 1 || significance <= 0 || significance >=1)
			stop("incorrect 'significance'")

		if (significance > 0.2) warning("'significance' is possibly incorrect")


		RESULT <- ifelse(alternative == "two.sided", (STATISTIC<qf(significance*0.5,2*ny,2*nx) || STATISTIC>qf(1-significance*0.5,2*ny,2*nx)),
		          ifelse(alternative == "greater",    STATISTIC>qf(1-significance,2*ny,2*nx),
		                                              STATISTIC<qf(significance,2*ny,2*nx)))

		RVAL <- list(statistic = STATISTIC, result = RESULT, alternative = nm_alternative,
			method = METHOD, data.name = DNAME)
		class(RVAL) <- "power.htest"
		return(RVAL)
	} else {

		PVAL <- ifelse(alternative == "two.sided", (0.5-abs(pf(STATISTIC, 2*ny, 2*nx)-0.5))*2,
		        ifelse(alternative == "greater",   1-pf(STATISTIC, 2*ny, 2*nx),
		                                             pf(STATISTIC, 2*ny, 2*nx)))


		RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative,
			method = METHOD, data.name = DNAME)
		class(RVAL) <- "htest"
		return(RVAL)
	}
}

