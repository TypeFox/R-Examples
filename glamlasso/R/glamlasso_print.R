#
#     Description of this R script:
#     R interface for glamlasso routines.
#
#     Intended for use with R.
#     Copyright (C) 2015 Adam Lund
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#

#' @title Print Function for objects of Class glamlasso
#'
#' @description This function will print some information about the glamlasso object.
#'  
#' @param x A glamlasso object
#' @param ... ignored
#' 
#' @examples  
#' \dontrun{
#' n1 <- 65; n2 <- 26; n3 <- 13; p1 <- 13; p2 <- 5; p3 <- 4
#' X1 <- matrix(rnorm(n1 * p1), n1, p1) 
#' X2 <- matrix(rnorm(n2 * p2), n2, p2) 
#' X3 <- matrix(rnorm(n3 * p3), n3, p3) 
#' Beta <- array(rnorm(p1 * p2 * p3) * rbinom(p1 * p2 * p3, 1, 0.1), c(p1 , p2, p3))
#' mu <- RH(X3, RH(X2, RH(X1, Beta)))
#' Y <- array(rnorm(n1 * n2 * n3, mu), dim = c(n1, n2, n3))
#' fit <- glamlasso(X1, X2, X3, Y, family = "gaussian", iwls = "exact")
#' fit
#' }
#' 
#' @method print glamlasso
#' @S3method print glamlasso
#' @author Adam Lund
#' @export
 
print.glamlasso <- function(x, ...) {
out <- data.frame(Df = x$df, lambda = x$lambda)
print(out)	
}
