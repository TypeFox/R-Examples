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

# Created: 2015.02.07
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

#' @title Higher order Cornish Fisher approximation.
#'
#' @description 
#'
#' Lee and Lin's Algorithm AS269 for higher order Cornish Fisher quantile approximation.
#'
#' @template details-cf
#'
#' @usage
#'
#' AS269(z,cumul,order.max=NULL,all.ords=FALSE)
#'
#' @param z the quantiles of the normal distribution. an atomic vector.
#' @param cumul the standardized cumulants of order 3, 4, ..., k. an atomic
#' vector.
#' @param order.max the maximum order approximation, must be greater than
#' \code{length(cumul)+2}.
#' We assume the cumulants have been adjusted to reflect that the random
#' variable has unit variance ('standardized cumulants')
#' @param all.ords a logical value. If \code{TRUE}, then results are returned
#' as a matrix, with a column for each order of the approximation. Otherwise
#' the results are a matrix with a single column of the highest order
#' approximation.
#' @return A matrix, which is, depending on \code{all.ords}, either with one column per 
#' order of the approximation, or a single column giving the maximum order
#' approximation. There is one row per value in \code{z}.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @note A warning will be thrown if any of the z are greater than 
#' 3.719017274 in absolute value; the traditional AS269 errors out in this
#' case.
#' @keywords distribution 
#' @seealso \code{\link{qapx_cf}}
#' @export 
#'
#' @template ref-AS269
#' @template ref-Jaschke
#'
#' @examples 
#' foo <- AS269(seq(-2,2,0.01),c(0,2,0,4))
#' # test with the normal distribution:
#' s.cumul <- c(0,0,0,0,0,0,0,0,0)
#' pv <- seq(0.001,0.999,0.001)
#' zv <- qnorm(pv)
#' apq <- AS269(zv,s.cumul,all.ords=FALSE)
#' err <- zv - apq
#'
#' # test with the exponential distribution
#' rate <- 0.7
#' n <- 18
#' # these are 'raw' cumulants'
#' cumul <- (rate ^ -(1:n)) * factorial(0:(n-1))
#' # standardize and chop
#' s.cumul <- cumul[3:length(cumul)] / (cumul[2]^((3:length(cumul))/2))
#' pv <- seq(0.001,0.999,0.001)
#' zv <- qnorm(pv)
#' apq <- cumul[1] + sqrt(cumul[2]) * AS269(zv,s.cumul,all.ords=TRUE)
#' truq <- qexp(pv, rate=rate)
#' err <- truq - apq
#' colSums(abs(err))
#'
#' # an example from Wikipedia page on CF, 
#' # \url{https://en.wikipedia.org/wiki/Cornish%E2%80%93Fisher_expansion}
#' s.cumul <- c(5,2)
#' apq <- 10 + sqrt(25) * AS269(qnorm(0.95),s.cumul,all.ords=TRUE)
#'
#' @template etc
AS269 <- function(z,cumul,order.max=NULL,all.ords=FALSE) {#FOLDUP
	if (is.null(order.max))
		order.max <- length(cumul) + 2
	if (order.max <= 2)
		return(z)
	stopifnot(order.max <= length(cumul)+2)
	if (any(abs(z) > 3.719017274))
		message('large z may result in inaccurate quantiles')

	# I use lower case letters for variables which are scalar
	# in z (or, rather, 'broadcast' across all z), and 
	# UPPER CASE for those which are the same size as 
	# z in one dimension or another. 
	# Otherwise, I have tried to make the variables match
	# the published AS269.
	nord <- order.max - 2

	# pre-line 10
	jidx <- 1:nord
	a <- ((-1)^jidx) * cumul / ((jidx+1) * (jidx+2))
	
	# line 10
	# prealloc H 
	n <- length(z)
	H <- matrix(0,nrow=n,ncol=3*nord)
	H[,1] <- - z
	H[,2] <- z * z - 1
	for (jj in 3:(dim(H)[2]))
		H[,jj] <- -(z * H[,jj-1] + (jj-1) * H[,jj-2])

	# line 20
	P <- matrix(0,nrow=n,ncol=3*nord*(nord+1)/2)

	# line 30
	D <- matrix(0,nrow=n,ncol=3*nord)
	DEL <- matrix(0,nrow=n,ncol=ifelse(all.ords,nord,1))

	D[,1] <- -a[1] * H[,2]
	DEL[,1] <- D[,1]
	P[,1] <- D[,1]
	P[,3] <- a[1]
	ja <- 0
	fac <- 1
	
	# main loop
	if (nord > 1) {
		for (jj in 2:nord) {#FOLDUP
			fac <- fac * jj
			ja <- ja + 3 * (jj-1)
			jb <- ja
			bc <- 1
			# calculate coefficients of Hermite polynomitals
			for (kk in 1:(jj-1)) {
				DD <- bc * D[,kk]
				aa <- bc * a[kk]
				jb <- jb - 3 * (jj - kk)
				for (ll in (1:(3*(jj-kk)))) {
					jbl <- jb + ll
					jal <- ja + ll
					P[,jal+1] <- P[,jal+1] + DD * P[,jbl]
					P[,jal+kk+2] <- P[,jal+kk+2] + aa * P[,jbl]
				}  # line 40
				bc <- bc * (jj - kk) / kk
			}  # line 50
			P[,ja + jj + 2] <- P[,ja + jj + 2] + a[jj]
			# calculate the adjustments
			D[,jj] <- 0
			for (ll in (2:(3*jj))) 
				D[,jj] <- D[,jj] - P[,ja + ll] * H[,ll-1]
			# line 60
			P[,ja+1] <- D[,jj]
			if (all.ords) 
				DEL[,jj] <- D[,jj] / fac
			else
				DEL <- DEL + (D[,jj] / fac)
		}  # line 70#UNFOLD
	}

	# the quantile approximations are then 
	# z + columns of DEL
	if (all.ords) 
		retval <- (z + t(apply(t(DEL),2,cumsum)))
	else
		retval <- z + DEL

	return(retval)
}#UNFOLD

#' @title Approximate quantile via Cornish-Fisher expansion.
#'
#' @description 
#'
#' Approximate the quantile function of a distribution via its cumulants.
#'
#' @details
#'
#' Given the cumulants of a probability distribution, we approximate the 
#' quantile function via a Cornish-Fisher expansion.
#'
#' @usage
#'
#' qapx_cf(p, raw.cumulants, support=c(-Inf,Inf), lower.tail = TRUE, log.p = FALSE)
#'
#' @param p where to evaluate the approximate distribution.
#' @param raw.cumulants an atomic array of the 1st through kth raw cumulants. The first 
#' value is the mean of the distribution, the second should
#' be the variance of the distribution, the remainder are raw cumulants.
#' @inheritParams dapx_gca
#' @return The approximate quantile at \code{p}.
#'
#' @keywords distribution 
#' @seealso \code{\link{dapx_gca}, \link{papx_gca}, \link{AS269}, \link{rapx_cf}}
#' @export 
#' @template ref-AS269
#' @template ref-Jaschke
#' @note 
#'
#' Monotonicity of the quantile function is not guaranteed.
#'
#' @examples 
#' # normal distribution:
#' pvals <- seq(0.001,0.999,length.out=501)
#' q1 <- qapx_cf(pvals, c(0,1,0,0,0,0,0))
#' q2 <- qnorm(pvals)
#' q1 - q2
#' @template etc
qapx_cf <- function(p,raw.cumulants,support=c(-Inf,Inf),lower.tail=TRUE,log.p=FALSE) {#FOLDUP
	order.max <- length(raw.cumulants)
	# this should be a standard routine:
	std.cumulants <- raw.cumulants / (raw.cumulants[2] ^ ((1:order.max)/2.0))
	z <- qnorm(p,lower.tail=lower.tail,log.p=log.p)
	if (order.max > 2) {
		gammas <- std.cumulants[3:order.max]
		z <- AS269(z=z,cumul=gammas,all.ords=FALSE)
	}

	# now convert back from standardized:
	retval <- raw.cumulants[1] + z * sqrt(raw.cumulants[2])
	# may have been converted to a matrix. fix that.
	dim(retval) <- dim(p)

	# support support
	if (is.finite(min(support))) {
		retval <- pmax(retval,min(support))
	}
	if (is.finite(max(support))) {
		retval <- pmin(retval,max(support))
	}

	return(retval)
}#UNFOLD

#' @title Approximate random generation via Cornish-Fisher expansion.
#'
#' @description 
#'
#' Approximate random generation via approximate quantile function.
#'
#' @details
#'
#' Given the cumulants of a probability distribution, we approximate the 
#' quantile function via a Cornish-Fisher expansion.
#'
#' @usage
#'
#' rapx_cf(n, raw.cumulants, support=c(-Inf,Inf))
#'
#' @param n number of observations. If 'length(n) > 1', the length is
#' taken to be the number required.
#' @inheritParams qapx_cf
#' @return A vector of approximate draws.
#'
#' @keywords distribution 
#' @seealso \code{\link{qapx_cf}}
#' @export 
#' @examples 
#' # normal distribution:
#' r1 <- rapx_cf(1000, c(0,1,0,0,0,0,0))
#' @template etc
rapx_cf <- function(n,raw.cumulants, support=c(-Inf,Inf)) {#FOLDUP
	p <- runif(n)
	retval <- qapx_cf(p=p,raw.cumulants=raw.cumulants,support=support)
	return(retval)
}#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
