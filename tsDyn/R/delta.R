## Copyright (C) 2005/2006  Antonio, Fabio Di Narzo
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.



#'delta test of conditional independence
#'
#'delta statistic of conditional independence and associated bootstrap test
#'
#'delta statistic of conditional independence and associated bootstrap test.
#'For details, see Manzan(2003).
#'
#'@aliases delta delta.test
#'@param x time series
#'@param m vector of embedding dimensions
#'@param d time delay
#'@param eps vector of length scales
#'@param B number of bootstrap replications
#'@return \code{delta} returns the computed delta statistic. \code{delta.test}
#'returns the bootstrap based 1-sided p-value.
#'@section Warning: Results are sensible to the choice of the window
#'\code{eps}. So, try the test for a grid of \code{m} and \code{eps} values.
#'Also, be aware of the course of dimensionality: m can't be too high for
#'relatively small time series. See references for further details.
#'@author Antonio, Fabio Di Narzo
#'@export
#'@seealso BDS marginal independence test: \code{\link[tseries]{bds.test}} in
#'package \pkg{tseries}
#'
#'Teraesvirta's neural network test for nonlinearity:
#'\code{\link[tseries]{terasvirta.test}} in package \pkg{tseries}
#'
#'delta test for nonlinearity: \code{\link{delta.lin.test}}
#'@references Sebastiano Manzan, Essays in Nonlinear Economic Dynamics, Thela
#'Thesis (2003)
#'@keywords ts
#'@examples
#'
#'delta(log10(lynx), m=3, eps=sd(log10(lynx)))
#'
delta <- function(x, m, d=1, eps) {
	if(m<2)
          stop("embedding dimension 'm' should be at least equal to 2")
	C <- d2(x, m=m+1, d=d, t=1, eps.min=eps)[1, 2:(m+2)]		#computed sample correlation integrals
	return( 1 - ( (C[m])^2 /  (C[m-1]*C[m+1]) ) )
}

#' @rdname delta
#' @export
delta.test <- function(x, m=2:3, d=1, eps = seq(0.5*sd(x),2*sd(x),length=4), B=49) {
	delta.b <- numeric(B)
	p.value <- matrix(NA,length(m),length(eps))
        for(j in 1:length(m)) for(i in 1:length(eps)){
		delta.c <- delta(x, m=m[j], d=d, eps=eps[i])
		for(b in 1:B)
                  delta.b[b] <- delta(sample(x), m=m[j], d=d, eps=eps[i])
		p.value[j,i] <- (1+sum(delta.b>=delta.c))/(1+B)
        }
	dimnames(p.value) <- list(m=m,eps=format(eps, digits=4))
	return(p.value)
}



#'delta test of linearity
#'
#'delta test of linearity based on conditional mutual information
#'
#'delta test of linearity based on conditional mutual information
#'
#'@aliases delta.lin delta.lin.test
#'@param x time series
#'@param m vector of embedding dimensions
#'@param d time delay
#'@param eps vector of length scales
#'@param B number of bootstrap replications
#'@return \code{delta.lin} returns the parametrically estimated delta statistic
#'for the given time series (assuming linearity). \code{delta.lin.test} returns
#'the bootstrap based 1-sided p-value. The test statistic is the difference
#'between the parametric and nonparametric delta estimators.
#'@author Antonio, Fabio Di Narzo
#'@references Sebastiano Manzan, Essays in Nonlinear Economic Dynamics, Thela
#'Thesis (2003)
#'@keywords ts
#'@export
#'@examples
#'
#'delta.lin(log10(lynx), m=3)
#'
delta.lin <- function(x, m, d=1) {
	V1 <- var(embedd(x, m=m+1, d=d))
	V2 <- var(embedd(x, m=m, d=d))
	tmp <- eigen(V1, symmetric=TRUE)$values[1] / eigen(V2,symmetric=TRUE)$values[1]
	return(1-tmp)
}

#' @rdname delta.lin
#' @export
delta.lin.test <- function(x, m=2:3, d=1, eps = seq(0.5*sd(x),2*sd(x),length=4), B=49) {
	mu <- function(x, m, eps)
          delta(x, m=m, d=d, eps=eps) - delta.lin(x, m=m, d=d)
	mu.c <- numeric()
	n <- length(x)
	ar.model <- ar(x)	#automatic AR order selection based on AIC
	mu.b <- numeric(B)
	p.value <- matrix(NA,length(m),length(eps))
        for(j in 1:length(m)) for(i in 1:length(eps)) {
          mu.c[i] <- mu(x, m[j], eps[i])
          for(b in 1:B) {
            xb <- arima.sim(n = n, list(ar=ar.model$ar),
                            rand.gen = function(n, ...) rnorm(n, 0, sqrt(ar.model$var.pred) ) )
            mu.b[b] <- mu(xb, m[j], eps[i])
          }
          p.value[j,i] <- (1+sum(mu.b>=mu.c[i]))/(1+B)
	}
	dimnames(p.value) <- list(m=m, eps=format(eps, digits=4))
	return(p.value)
}
