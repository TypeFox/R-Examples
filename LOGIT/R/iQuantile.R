# Hilbe, J.M., Practical Guide to Logistic Regression 2015
# Rafael de Souza, Eotvos Lorand Univ.
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as published by
#the Free Software Foundation.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
#' @title  iQuantile
#' @description Produces indices in x  of percentile steps by 1/breaks and respective percentiles.
#' @aliases iQuantile
#' @usage iQuantile(x, breaks=15)
#' @importFrom stats quantile sd
#' @format \describe{
#' \item{x}{
#' The function has two arguments: a  list and numnber of breaks.}
#' }
#' @details iQuantile is used internally by the  hlGOF.test.
#' @param x list
#' @param breaks number of breaks
#' @return numeric
#'@examples
#'library(LOGIT)
#'mod <- rnorm(100,0,1)
#'iQuantile(mod)
#' @author Joseph M. Hilbe, Arizona State University.
#'
#' @references Hilbe, J. M. (2015), Practical Guide to Logistic Regression, Chapman & Hall/CRC.
#'
#' Hilbe, J. M. (2009), Logistic Regression Models, Chapman & Hall/CRC.
#' @keywords models
#' @export
#'
iQuantile <- function (x, breaks=15) {
  #indices in x[] of percentile steps by 1/breaks
  xo <- order(x)  #sort indices
  n <- length(x)
  r <- rep(0, breaks+1)
  r[1]<- 1
  r[breaks+1]<- n
  r[2:breaks]<- round(1:(breaks-1)*n/breaks)
  return(list(index=r,cuts=x[xo[r]]))
}
