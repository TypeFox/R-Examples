# Copyright 2016-2016 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

# This file is part of BWStest.
#
# BWStest is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BWStest is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BWStest.  If not, see <http://www.gnu.org/licenses/>.

#' @title Perform the Baumgartner-Weiss-Schindler hypothesis test.
#'
#' @description 
#'
#' Perform the Baumgartner-Weiss-Schindler hypothesis test.
#'
#' @param x a vector of the first sample.
#' @param y a vector of the first sample.
#' @param method a character string specifying the test statistic to use.
#' should be one of the following:
#' \describe{
#' \item{default}{This is \dQuote{Hobson's choice}, which uses the classical
#' BWS test for two-sided alternative, but Neuhauser for one sided
#' alternatives.}
#' \item{BWS}{Use the classical BWS test.}
#' \item{Neuhauser}{Use Neuhauser's test.}
#' \item{B1}{Use Murakami's \eqn{B_1}{B1} test.}
#' \item{B2}{Use Murakami's \eqn{B_2}{B2} test, which is exactly Neuhauser's test.}
#' \item{B3}{Use Murakami's \eqn{B_3}{B3} test.}
#' \item{B4}{Use Murakami's \eqn{B_4}{B4} test.}
#' \item{B5}{Use Murakami's \eqn{B_5}{B5} test.}
#' }
#' Only Neuhauser's test supports one-sided alternatives.
#' @param alternative a character string specifying the alternative hypothesis,
#'        must be one of \dQuote{two.sided} (default), \dQuote{greater} or
#'        \dQuote{less}. You can specify just the initial letter.
#'        \dQuote{greater} corresponds to testing whether the survival function
#'        of \code{x} is greater than that of \code{y}; equivalently one can
#'        think of this as \code{x} being \sQuote{greater} than \code{y}
#'        in the sense of first order stochastic dominance.
#' @return Object of class \code{htest}, a list of the test statistic,
#' the p-value, and the \code{method} noted.
#' @keywords htest
#' @seealso \code{\link{bws_test}}, \code{\link{bws_stat}},
#' \code{\link{murakami_stat}}, \code{\link{murakami_cdf}}.
#' @template etc
#' @template ref-bws
#' @note The code will happily compute Murakami's \eqn{B_3} through \eqn{B_5}
#' for large sample sizes, even though nominal coverage is \emph{not} achieved.
#' A warning will be thrown. User assumes all risk relying on results from this
#' function.
#' @examples 
#'
#' # under the null
#' set.seed(123)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' hval <- bws_test(x,y)
#' 
#' # under the alternative
#' set.seed(123)
#' x <- rnorm(100)
#' y <- rnorm(100,mean=1.0)
#' hval <- bws_test(x,y)
#' show(hval)
#' stopifnot(hval$p.value < 0.05)
#' 
#' # under the alternative with a one sided test.
#' set.seed(123)
#' x <- rnorm(100)
#' y <- rnorm(100,mean=0.7)
#' hval <- bws_test(x,y,alternative='less')
#' show(hval)
#' stopifnot(hval$p.value < 0.01)
#'
#' hval <- bws_test(x,y,alternative='greater')
#' stopifnot(hval$p.value > 0.99)
#'
#' hval <- bws_test(x,y,alternative='two.sided')
#' stopifnot(hval$p.value < 0.05)
#'
#' @rdname bws_test
#' @export
bws_test <- function(x,y,
										 method=c('default','BWS','Neuhauser','B1','B2','B3','B4','B5'),
										 alternative=c("two.sided","greater","less")) 
{
	method <- match.arg(method)
	alternative <- match.arg(alternative)
	dname <- paste(deparse(substitute(x)),'vs.',deparse(substitute(y)))
	x <- x[!is.na(x)]
	y <- y[!is.na(y)]
	nx <- length(x)
	ny <- length(y)
	if (max(nx,ny) <= 8) {
		warning('A permutation test would likely make more sense.')
	} else if (min(nx,ny) <= 10) {
		warning('Small, imbalanced, sample size may cause loss of nominal coverage.')
	}

	if (method == 'Neuhauser') { method <- 'B2' }  # for simplicity

	if (alternative == 'two.sided') {
		if (method == 'default') { method <- 'BWS' }  
		switch(method,
					 BWS={
							method <- "two-sample BWS test"
							bval <- bws_stat(x,y)
							names(bval) <- "B"
							pval <- bws_cdf(bval,lower_tail=FALSE)
					},
					{  # default is parametrized by the number.
						flavor <- as.numeric(gsub('^B(\\d)$','\\1',method))
					 	method <- "two-sample Murakami test"
					 	bval <- murakami_stat(x,y,flavor=flavor)
					 	names(bval) <- sprintf('B_%d',flavor)
					 	pval <- murakami_cdf(bval,n1=nx,n2=ny,flavor=flavor,lower_tail=FALSE)
					})
	} else {
		if (method == 'default') { method <- 'B2' }  
		stopifnot(method %in% c('B2'))
		method <- "two-sample Neuhauser/Murakami test"
		switch(alternative,
					 greater={
						 bval <- murakami_stat(x,y,flavor=2L)
						 pval <- murakami_cdf(bval,n1=nx,n2=ny,flavor=2L)
					 },
					 less={
						 bval <- murakami_stat(y,x,flavor=2L)
						 pval <- murakami_cdf(bval,n1=ny,n2=nx,flavor=2L)
					 })
		names(bval) <- "B_2"
	}
	zeta <- 0
	names(zeta) <- "difference in survival functions"
	if ((method %in% c('B3','B4','B5')) && (max(nx,ny) > 12)) {
		warning('Nominal coverage for B3 through B5 is *not* achieved for larger sample sizes. Use at your own risk!')
	}

	retval <- list(statistic = bval, 
								 p.value = pval,
								 alternative = alternative,
								 null.value = zeta,
								 method = method, 
								 data.name = dname)
	class(retval) <- "htest"
	return(retval)
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
