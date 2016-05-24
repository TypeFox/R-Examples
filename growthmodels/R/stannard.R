##
##  Stannard growth model
##
##  Created by Daniel Rodríguez Pérez on 28/8/2013.
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

#' Stannard growth model
#'
#' Computes the Stannard growth model
#' \deqn{ y(t) = \alpha \left[ 1 + exp(-(\beta + k t)/m) \right]^{-m}}{ y(t) = \alpha *( 1 + exp(-(beta + k * t)/m))^(-m) }
#' 
#' @param t time
#' @param x size
#' @param alpha upper asymptote
#' @param beta growth displacement
#' @param k growth rate 
#' @param m slope of growth 
#' 
#' @usage stannard(t, alpha, beta, k, m)
#' 
#' @examples
#' growth <- stannard(0:10, 1, .2, .1, .5)
#' 
#' @references
#' A. Khamiz, Z. Ismail, and A. T. Muhammad, "Nonlinear growth models for
#' modeling oil palm yield growth," Journal of Mathematics and Statistics,
#' vol. 1, no. 3, p. 225, 2005.
#' 
#' @rdname stannard
#' @export stannard
#' @aliases stannard
stannard <- function(t, alpha, beta, k, m) {
  result <- alpha * ( 1 + exp(-(beta + k * t)/m) )^(-m)
  return(result)
}

#' @examples
#' # Calculate inverse function
#' time <- stannard.inverse(growth, 1, .2, .1, .5)
#' 
#' @rdname stannard
#' @export stannard.inverse
#' @aliases stannard.inverse
stannard.inverse <- function(x, alpha, beta, k, m) {
  result <- - (beta + m * log((alpha / x)^(1 / m) - 1)) / k
  return(result)
}
