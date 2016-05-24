##
##  Chapman-Richards growth model
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

#' Chapman-Richards growth model
#'
#' Computes the Chapman-Richards growth model and its inverse
#' \deqn{ y(t) = \alpha (1 - \beta exp(-k t)^{1/(1-m)}) }{ y(t) = \alpha * (1 - \beta * exp(-k * t)^{1/(1-m)}) }
#' 
#' @param t time
#' @param x  size
#' @param alpha upper asymptote
#' @param beta growth range 
#' @param k growth rate
#' @param m slope of growth 
#' 
#' @usage chapmanRichards(t, alpha, beta, k, m)
#' 
#' @examples
#' growth <- chapmanRichards(0:10, 10, 0.5, 0.3, 0.5)
#' 
#' @references
#' D. Fekedulegn, M. Mac Siurtain, and J. Colbert, "Parameter estimation of
#' nonlinear growth models in forestry," Silva Fennica, vol. 33, no. 4, pp.
#' 327-336, 1999.
#' 
#' @rdname chapmanRichards
#' @export chapmanRichards
#' @aliases chapmanRichards
chapmanRichards <- function(t, alpha, beta, k, m) {
  result <- alpha * (1 - beta * exp(-k * t))^(1 / (1 - m))
  return(result)
}

#' @examples
#' # Calculate inverse function
#' time <- chapmanRichards.inverse(growth, 10, 0.5, 0.3, 0.5)
#' 
#' @rdname chapmanRichards
#' @export chapmanRichards.inverse
#' @aliases chapmanRichards.inverse
chapmanRichards.inverse <- function(x, alpha, beta, k, m) {
  result <- - log((1 - (x / alpha)^(1-m)) / beta) / k
  return(result)
}
