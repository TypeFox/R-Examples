##
##  Weibull growth model
##
##  Created by Daniel Rodríguez Pérez on 28/7/2013.
##
##  Copyright (c) 2013 Daniel Rodríguez Pérez.
##
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>
##

#' Weibull growth model
#'
#' Computes the Weibull growth model
#' \deqn{ y(t) = \alpha - \beta exp(-k * t^m) }{ y(t) = \alpha - \beta * exp(-k * t^m) }
#' 
#' @param t time
#' @param x size
#' @param alpha upper asymptote
#' @param beta growth range 
#' @param k growth rate
#' @param m slope of growth 
#' 
#' @usage weibull(t, alpha, beta, k, m)
#' 
#' @examples
#' growth <- weibull(0:10, 10, 0.5, 0.3, 0.5)
#' 
#' @references
#' D. Fekedulegn, M. Mac Siurtain, and J. Colbert, "Parameter estimation of
#' nonlinear growth models in forestry," Silva Fennica, vol. 33, no. 4, pp.
#' 327-336, 1999.
#' 
#' @rdname weibull
#' @export weibull
#' @aliases weibull
weibull <- function(t, alpha, beta, k, m) {
  result <- alpha - beta * exp(-k * t^m);
  return(result)
}

#' @examples
#' # Calculate inverse function
#' time <- weibull.inverse(growth, 10, 0.5, 0.3, 0.5)
#' 
#' @rdname weibull
#' @export weibull.inverse
#' @aliases weibull.inverse
weibull.inverse <- function(x, alpha, beta, k, m) {
  result <- ((-1/k) * log((alpha - x)/beta))^(1/m)
  return(result)
}
