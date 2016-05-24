# Copyright 2015-2015 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of PDQutils.
#
# PDQutils is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PDQutils is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with PDQutils.  If not, see <http://www.gnu.org/licenses/>.

# Created: 2015.02.08
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# utilities or dealing with moments and cumulants

#' @title Convert moments to raw cumulants.
#'
#' @description 
#'
#' Conversion of a vector of moments to raw cumulants.
#'
#' @details
#'
#' The 'raw' cumulants \eqn{\kappa_i}{kappa_i} are connected
#' to the 'raw' (uncentered) moments, \eqn{\mu_i'}{mu'_i} via
#' the equation
#' \deqn{\kappa_n = \mu_n' - \sum_{m=1}^{n-1} {n-1 \choose m-1} \kappa_m \mu_{n-m}'}
#'
#' Note that this formula also works for central moments, assuming
#' the distribution has been normalized to zero mean.
#'
#' @usage
#'
#' moment2cumulant(moms)
#'
#' @param moms a vector of the moments. The first element is the first moment.
#' If centered moments are given, the first cumulant shall be zero.
#' @return a vector of the cumulants.
#'
#' @keywords distribution 
#' @seealso \code{\link{cumulant2moment}}
#' @export 
#'
#' @examples 
#' # normal distribution, mean 0, variance 1
#' n.cum <- moment2cumulant(c(0,1,0,3,0,15))
#' # normal distribution, mean 1, variance 1
#' n.cum <- moment2cumulant(c(1,2,4,10,26))
#' # exponential distribution
#' lambda <- 0.7
#' n <- 1:6
#' e.cum <- moment2cumulant(factorial(n) / (lambda^n))
#' @template etc
moment2cumulant <- function(moms) {
	kappa <- moms
	if (length(kappa) > 1) {
		for (nnn in 2:length(kappa)) {
			mmm <- 1:(nnn-1)
			kappa[nnn] <- moms[nnn] - sum(choose(nnn-1,mmm-1) * kappa[mmm] * moms[nnn-mmm])
		}
	}
	return(kappa)
}

#' @title Convert raw cumulants to moments.
#'
#' @description 
#'
#' Conversion of a vector of raw cumulatnts to moments.
#'
#' @details
#'
#' The 'raw' cumulants \eqn{\kappa_i}{kappa_i} are connected
#' to the 'raw' (uncentered) moments, \eqn{\mu_i'}{mu'_i} via
#' the equation
#' \deqn{\mu_n' = \kappa_n + \sum_{m=1}^{n-1} {n-1 \choose m-1} \kappa_m \mu_{n-m}'}
#'
#' @usage
#'
#' cumulant2moment(kappa)
#'
#' @param kappa a vector of the raw cumulants. The first element is the first cumulant,
#' which is also the first moment.
#' @return a vector of the raw moments.
#'
#' @keywords distribution 
#' @seealso \code{\link{moment2cumulant}}
#' @export 
#'
#' @examples 
#' # normal distribution, mean 0, variance 1
#' n.mom <- cumulant2moment(c(0,1,0,0,0,0))
#' # normal distribution, mean 1, variance 1
#' n.mom <- cumulant2moment(c(1,1,0,0,0,0))
#' @template etc
cumulant2moment <- function(kappa) {
	moms <- kappa
	if (length(moms) > 1) {
		for (nnn in 2:length(kappa)) {
			mmm <- 1:(nnn-1)
			moms[nnn] <- kappa[nnn] + sum(choose(nnn-1,mmm-1) * kappa[mmm] * moms[nnn-mmm])
		}
	}
	return(moms)
}

# central moments to standardized moments
# assumes given _zeroth_ through kth centered moment.
# divides through by the standard deviation to the
# correct power.
# not needed for now..
#central2std <- function(mu.cent) {
	#mu.std <- mu.cent / (mu.cent[3] ^ ((0:(length(mu.cent)-1))/2))
	#return(mu.std)
#}


#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
