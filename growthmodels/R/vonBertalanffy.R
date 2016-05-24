##
##  von Bertalanffy growth model
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

#' von Bertalanffy growth model
#'
#' Computes the von Bertalanffy growth model
#' \deqn{ y(t) = (\alpha^(1-m) - \beta * exp(-k t))^(1/(1-m)) }{ y(t) = (\alpha^(1-m) - \beta * exp(-k * t))^(1/(1-m)) }
#' 
#' @param t time
#' @param x size
#' @param alpha upper asymptote
#' @param beta growth range 
#' @param k growth rate
#' @param m slope of growth 
#' 
#' @usage vonBertalanffy(t, alpha, beta, k, m)
#' 
#' @examples
#' growth <- vonBertalanffy(0:10, 10, 0.5, 0.3, 0.5)
#' 
#' @references
#' D. Fekedulegn, M. Mac Siurtain, and J. Colbert, "Parameter estimation of
#' nonlinear growth models in forestry," Silva Fennica, vol. 33, no. 4, pp.
#' 327-336, 1999.
#' 
#' @rdname vonBertalanffy
#' @export vonBertalanffy
#' @aliases vonBertalanffy
vonBertalanffy <- function(t, alpha, beta, k, m) {
  result <- alpha^(1 - m) - beta * exp(-k * t)
  result <- result^(1 / (1 - m))
  return(result)
}

#' @examples
#' # Calculate inverse function
#' time <- vonBertalanffy.inverse(growth, 10, 0.5, 0.3, 0.5)
#' 
#' @rdname vonBertalanffy
#' @export vonBertalanffy.inverse
#' @aliases vonBertalanffy.inverse
vonBertalanffy.inverse <- function(x, alpha, beta, k, m) {
  result <- -log((alpha^(1 - m) - x^(1 - m))/beta) / k
  return(result)
}
