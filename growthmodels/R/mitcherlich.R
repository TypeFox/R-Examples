##
##  Mitcherlich exponential growth model
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

#' Mitcherlich growth model
#'
#' Computes the Mitcherlich growth model
#' \deqn{ y(t) = (\alpha - \beta k^t)}{ y(t) = \alpha - \beta * k^t}
#' 
#' @param t time
#' @param x size
#' @param alpha upper asymptote
#' @param beta growth range 
#' @param k growth rate 
#' 
#' @usage mitcherlich(t, alpha, beta, k)
#' 
#' @examples
#' growth <- mitcherlich(0:10, 10, 0.5, 0.3)
#' 
#' @references
#' D. Fekedulegn, M. Mac Siurtain, and J. Colbert, "Parameter estimation of
#' nonlinear growth models in forestry," Silva Fennica, vol. 33, no. 4, pp.
#' 327-336, 1999.
#' 
#' @rdname mitcherlich
#' @export mitcherlich
#' @aliases mitcherlich
mitcherlich <- function(t, alpha, beta, k) {
  result <- alpha - beta * k^t
  return(result)
}

#' @examples
#' # Calculate inverse function
#' time <- mitcherlich.inverse(growth, 10, 0.5, 0.3)
#' 
#' @rdname mitcherlich
#' @export mitcherlich.inverse
#' @aliases mitcherlich.inverse
mitcherlich.inverse <- function(x, alpha, beta, k) {
  result <- log((alpha - x) / beta) / log(k)
  return(result)
}
