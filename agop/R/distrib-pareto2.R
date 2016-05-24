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


#' @title
#' Pareto Type-II (Lomax) Distribution
#' 
#' @description
#' Density, cumulative distribution function,
#' quantile function, and random generation for the 
#' Pareto Type-II (Lomax)  
#' distribution with shape
#' parameter \eqn{k>0} and scale parameter \eqn{s>0}.
#' 
#' @details
#' If \eqn{X\sim\mathrm{P2}(k,s)}{X~P2(k,s)},
#' then \eqn{\mathrm{supp}\,X=[0,\infty)}{supp X=[0,\infty)}.
#' The c.d.f. for \eqn{x\ge 0} is given by \deqn{F(x)=1-s^k/(s+x)^k}
#' and the density by \deqn{f(x)=k s^k/(s+x)^{k+1}.}
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n integer; number of observations
#' @param k vector of shape parameters, \eqn{k>0}
#' @param s vector of scale parameters, \eqn{s>0}
#' @param lower.tail logical; if \code{TRUE} (default),
#' probabilities are \eqn{P(X \le x)}, and \eqn{P(X > x)} otherwise
#' @return 
#' numeric vector;
#' \code{dpareto2} gives the density,
#' \code{ppareto2} gives the cumulative distribution function,
#' \code{qpareto2} calculates the quantile function,
#' and \code{rpareto2} generates random deviates.
#' 
#' @export
#' @rdname Pareto2
#' @family distributions
#' @family Pareto2
rpareto2 <- function(n, k=1, s=1)
{
   stopifnot(is.numeric(k), k > 0)
   stopifnot(is.numeric(s), s > 0)
   # n checked by runif
	s*((runif(n)^(-1/k)) - 1)
}



#' @export
#' @rdname Pareto2
ppareto2 <- function(q, k=1, s=1, lower.tail=TRUE)
{
   stopifnot(is.numeric(k), k > 0)
   stopifnot(is.numeric(s), s > 0)
   stopifnot(is.numeric(q))
	ret <- ifelse(q<0, 0, (1-(s/(s+q))^k))
   if (identical(lower.tail[1], FALSE))
      ret <- 1-ret
   else
      ret
}



#' @export
#' @rdname Pareto2
qpareto2 <- function(p, k=1, s=1, lower.tail=TRUE)
{
   stopifnot(is.numeric(k), k > 0)
   stopifnot(is.numeric(s), s > 0)
   stopifnot(is.numeric(p))
   if (identical(lower.tail[1], FALSE))
      p <- 1-p
	ifelse(p<0 | p>1, NaN, s*((1-p)^(-1/k)-1))
}



#' @export
#' @rdname Pareto2
dpareto2 <- function(x, k=1, s=1)
{
   stopifnot(is.numeric(k), k > 0)
   stopifnot(is.numeric(s), s > 0)
   stopifnot(is.numeric(x))
	ifelse(x<=0, 0, k/(s+x)*(s/(s+x))^k)
}
