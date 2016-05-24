##
##  Blumberg growth model
##
##  Created by Daniel Rodríguez Pérez on 14/9/2013.
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

#' Blumberg growth model
#'
#' Computes the Blumberg growth model and its inverse
#' \deqn{ y(t) = \frac{\alpha * (t + t_0)^m}{w_0 + (t + t_0)^m}}{y(t) = (\alpha * (t - t_0)^m)/(w_0 + (t - t_0)^m)}
#' 
#' @param t time
#' @param x size
#' @param alpha upper asymptote
#' @param w0 a reference value at t = t0
#' @param m slope of growth 
#' @param t0 time shift (default 0)
#' 
#' @examples
#' growth <- blumberg(0:10, 10, 2, 0.5)
#' 
#' @references
#' A. Tsoularis and J. Wallace, "Analysis of logistic growth models.,"
#' Math Biosci, vol. 179, no. 1, pp. 21-55, Jul. 2002.
#' 
#' @rdname blumberg
#' @export blumberg
#' @aliases blumberg
blumberg <- function(t, alpha, w0, m, t0 = 0) {
  result <- (alpha * (t + t0)^m) / (w0 + (t + t0)^m)
  return(result)
}

#' @examples
#' # Calculate inverse function
#' time <- blumberg.inverse(growth, 12, 2, 0.5) 
#' 
#' @rdname blumberg
#' @export blumberg.inverse
#' @aliases blumberg.inverse
blumberg.inverse <- function(x, alpha, w0, m, t0 = 0) {
  result <- (x * w0 / (alpha - x))^(1/m) - t0
  return(result)
}
