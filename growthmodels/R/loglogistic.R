##
##  Log-logistic growth model
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

#' Log-logistic growth model
#'
#' Computes the Log-logistic growth model
#' \deqn{ y(t) = \frac{\alpha}{1 + \beta exp(-k log(t)}}{ y(t) = \alpha/(1 + \beta * exp(-k * log(t))}
#' 
#' @param t time
#' @param x size
#' @param alpha upper asymptote
#' @param beta growth range 
#' @param k growth rate 
#' 
#' @usage loglogistic(t, alpha, beta, k)
#' 
#' @examples
#' growth <- loglogistic(0:10, 10, 0.5, 0.3)
#' 
#' @references
#' A. Khamiz, Z. Ismail, and A. T. Muhammad, "Nonlinear growth models for
#' modeling oil palm yield growth," Journal of Mathematics and Statistics,
#' vol. 1, no. 3, p. 225, 2005.
#' 
#' @rdname loglogistic
#' @export loglogistic
#' @aliases loglogistic
loglogistic <- function(t, alpha, beta, k) {
  t[t < 0] <- NaN
  result   <- logistic(log(t), alpha, beta, k)
  return(result)
}

#' @examples
#' # Calculate inverse function
#' time <- loglogistic.inverse(growth, 10, 0.5, 0.3)
#' 
#' @rdname loglogistic
#' @export loglogistic.inverse
#' @aliases loglogistic.inverse
loglogistic.inverse <- function(x, alpha, beta, k) {
  result <- exp(logistic.inverse(x, alpha, beta, k))
  return(result)
}

