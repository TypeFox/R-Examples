##
##  Schnute growth model
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

#' Schnute growth model
#'
#' Computes the Schnute growth model
#' \deqn{ y(t) =  \left[ r_0 + \beta exp(k t) \right]^m }{ y(t) = (r_0 + \beta * exp(k * t))^m }
#' 
#' @param t time
#' @param x size
#' @param r0 reference value
#' @param beta growth displacement
#' @param k growth rate 
#' @param m slope of growth 
#' 
#' @usage schnute(t, r0, beta, k, m)
#' 
#' @examples
#' growth <- schnute(0:10, 10, 5, .5, .5)
#' 
#' @references
#' A. Khamiz, Z. Ismail, and A. T. Muhammad, "Nonlinear growth models for
#' modeling oil palm yield growth," Journal of Mathematics and Statistics,
#' vol. 1, no. 3, p. 225, 2005.
#' 
#' @rdname schnute
#' @export schnute
#' @aliases schnute
schnute <- function(t, r0, beta, k, m) {
  result <- (r0 + beta * exp(k * t))^m
  return(result)
}

#' @examples
#' # Calculate inverse function
#' time <- schnute.inverse(growth, 10, 5, .5, .5)
#' 
#' @rdname schnute
#' @export schnute.inverse
#' @aliases schnute.inverse
schnute.inverse <- function(x, r0, beta, k, m) {
  result <- log((x^(1 / m) - r0) / beta) / k
  return(result)
}
