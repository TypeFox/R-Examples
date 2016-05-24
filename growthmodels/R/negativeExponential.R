##
##  Negative exponential growth model
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

#' Negative exponential growth model
#'
#' Computes the negative exponential growth model
#' \deqn{ y(t) = \alpha ( 1 - exp(-k t))}{ y(t) = \alpha * ( 1 - exp(-k * t))}
#' 
#' @param t time
#' @param x size
#' @param alpha upper asymptote
#' @param k growth rate 
#' 
#' @usage negativeExponential(t, alpha, k)
#' 
#' @examples
#' growth <- negativeExponential(0:10, 1, 0.3)
#' 
#' @references
#' D. Fekedulegn, M. Mac Siurtain, and J. Colbert, "Parameter estimation of
#' nonlinear growth models in forestry," Silva Fennica, vol. 33, no. 4, pp.
#' 327-336, 1999.
#' 
#' @rdname negativeExponential
#' @export negativeExponential
#' @aliases negativeExponential
negativeExponential <- function(t, alpha, k) {
  result <- alpha * (1.0 - exp(-k * t))
  return(result)
}

#' @examples
#' # Calculate inverse function
#' time <- negativeExponential.inverse(growth, 10, 0.3)
#' 
#' @rdname negativeExponential
#' @export negativeExponential.inverse
#' @aliases negativeExponential.inverse
negativeExponential.inverse <- function(x, alpha, k) {
  result <- -log(1 - x/alpha) / k
  return(result)
}
