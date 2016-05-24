##
##  Gompertz exponential growth model
##
##  Created by Daniel Rodríguez Pérez on 27/7/2013.
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

#' Gompertz growth model
#'
#' Computes the Gompertz growth model and its inverse
#' \deqn{ y(t) = \alpha exp(-\beta exp(-k^t))}{ y(t) = \alpha * exp(-\beta * exp(-k^t))}
#' 
#' @param t time
#' @param x size
#' @param alpha upper asymptote
#' @param beta growth displacement
#' @param k growth rate 
#' 
#' @usage gompertz(t, alpha, beta, k)
#' 
#' @examples
#' growth <- gompertz(0:10, 10, 0.5, 0.3)
#' 
#' @references
#' D. Fekedulegn, M. Mac Siurtain, and J. Colbert, "Parameter estimation of
#' nonlinear growth models in forestry," Silva Fennica, vol. 33, no. 4, pp.
#' 327-336, 1999.
#' 
#' @rdname gompertz
#' @export gompertz
#' @aliases gompertz
gompertz <- function(t, alpha, beta, k) {
  result <- alpha * exp(-beta * exp(-k * t));
  return(result)
}

#' @examples
#' # Calculate inverse function
#' time <- gompertz.inverse(growth, 10, 0.5, 0.3)
#' 
#' @rdname gompertz
#' @export gompertz.inverse
#' @aliases gompertz.inverse
gompertz.inverse <- function(x, alpha, beta, k) {
  result <- - log(-log(x / alpha) / beta) / k
  return(result)
}
