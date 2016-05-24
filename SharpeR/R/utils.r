# Copyright 2012-2014 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

# This file is part of SharpeR.
#
# SharpeR is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SharpeR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SharpeR.  If not, see <http://www.gnu.org/licenses/>.

# env var:
# nb: 
# see also:
# todo:
# changelog: 
#
# Created: 2012.05.19
# Copyright: Steven E. Pav, 2012-2013
# Author: Steven E. Pav
# Comments: Steven E. Pav

# matrix utils#FOLDUP
# asymetric toeplitz, like Matlab's; 
# allows you to give the first row vector and the first column vector.
.atoeplitz <- function(rv,cv=rv) {
	if (!is.vector(rv)) 
		stop("'rv' is not a vector")
	if (!is.vector(cv)) 
		stop("'cv' is not a vector")
	nc <- length(rv)
	nr <- length(cv)
	if (rv[1] != cv[1])
		stop("'rv' and 'cv' must match at first element")

	A <- matrix(raw(), nr, nc)
	retval <- matrix(ifelse((col(A) > row(A)),rv[abs(col(A) - row(A)) + 1L],cv[abs(col(A) - row(A)) + 1L]), nr, nc)
	return(retval)
}

# the mean gram 
#.mean.gram <- function(X) { (t(X) %*% X) / dim(X)[1] }

# quadratic forms
.qform <- function(bread,meat) { # nocov start
	return(t(bread) %*% (meat %*% bread))
} # nocov end
# and 'outer' 
.qoform <- function(bread,meat) {
	return(bread %*% (meat %*% t(bread)))
}
#UNFOLD

# stats utils#FOLDUP
# convert p-value for a 1-sided test into one for a two-sided test.
.oneside2two <- function(pv) {
	retval <- 1 - 2 * abs(0.5 - pv)
	return(retval)
}
#UNFOLD

# class utils#FOLDUP
# infer the total time from an xts object,
# as a difftime
.infer_delt_xts <- function(anxts) {
	TEO <- time(anxts)
	delt <- difftime(TEO[length(TEO)],TEO[1],units='days')
	return(delt)
}
# infer the observations per epoch from an xts object
.infer_ope_xts <- function(anxts) {
	delt <- .infer_delt_xts(anxts)
	n.row <- dim(anxts)[1]
	days.per.row <- as.numeric(delt) / (n.row - 1)
	ope <- 365.25 / days.per.row
	return(ope)
}
#UNFOLD

# annualize and deannualize a Sharpe Ratio#FOLDUP
# ' @param sr the Sharpe Ratio, in per sqrt(epoch) units.
# ' @param ope the number of observations per epoch. no default here.
# ' @return the annualized Sharpe ratio, in per sqrt(epoch) units.
# ' @aliases .deannualize
# ' @export
.annualize <- function(sr, ope) {
  return(sr * sqrt(ope))  
}
# ' @export
.deannualize <- function(sr.pa, ope) {
  return(sr.pa / sqrt(ope))
}

# for T^2: 
.annualize2 <- function(T2, ope) { # nocov start
  return(T2 * ope)  
} # nocov end
.deannualize2 <- function(T2.pa, ope) { # nocov start
  return(T2.pa / ope)
} # nocov end
#UNFOLD

# conversions
# converting T2 <-> F#FOLDUP
# convert hotelling T2 to F statistic
.T2_to_F <- function(T2, p, n) {
	return(T2 * (n - p) / (p * (n - 1)))
}
# derivative of same
.d_T2_to_F <- function(T2, p, n) {
	return((n - p) / (p * (n - 1)))
}
# don't need this?
.logT2_to_logF <- function(logT2, p, n) {
	return(logT2 + log(.T2_to_F(1, p, n)))  # correct but mildly inefficient
}
# convert F statistic to hotelling T2
.F_to_T2 <- function(F, p, n) {
	return(F * (p * (n - 1)) / (n - p))
}
#.logF_to_logT2 <- function(logF, p, n) {
	#return(logF + .F_to_T2(1, p, n))  # correct but mildly inefficient
#}
#UNFOLD

# converting T2 <-> sropt#FOLDUP
# convert hotelling T2 to maximal SR statistic, SR^*
.T2_to_sropt <- function(T2, n) {
	return(sqrt(T2 / n))
}
# derivative of same
.d_T2_to_sropt <- function(T2, n) {
	return(0.5 / sqrt(T2 * n))
}
.logT2_to_logsropt <- function(logT2, n) {
	return(0.5 * (logT2 - log(n)))
}
# convert sropt to T2; 
# 2FIX: this is a total cluster fuck b/c negative sropt
# are 'out of range'; silently mark them up to zero here
# and pretend nothing happened? crap.
.sropt_to_T2 <- function(sropt, n) {
	return(n * pmax(0,sropt)^2)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
