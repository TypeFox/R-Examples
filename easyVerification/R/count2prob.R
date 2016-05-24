# count2prob Convert Ensemble Counts to Probabilities
#
#     Copyright (C) 2016 MeteoSwiss
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
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#' Convert Ensemble Counts to Probabilities
#' 
#' Using plotting positions as described in Wilks (2011), counts of 
#' occurences per forecast category are converted to probabilities
#' of occurence. For ensembles of size 1 (e.g. verifying observations), 
#' the count vector is returned unaltered (corresponding to occurence 
#' probabilities of 0 or 1).
#' 
#' 
#' @param x input matrix of counts from \code{\link{convert2prob}}
#' @param type selection of plotting positions (default to 3, see Types)
#' 
#' @section Types:
#' 
#' The types characterize the plotting positions as specified in Wilks (2011). The 
#' plotting positions are computed using the following relationship:
#' 
#' \deqn{p(x_i) = \frac{i + 1 - a}{n + 1 - a}}{p(x_i) = (i + 1 - a)/(n + 1 - a)}
#' 
#' where i is the number of ensemble members not exceeding x, and n is the
#' number of ensemble members. The types are characterized as follows:
#' 
#' \tabular{clc}{
#'   type \tab description \tab a \cr
#'   1 \tab Weibull \tab 0 \cr
#'   2 \tab Bernard and Bos-Levenbach \tab 0.3 \cr
#'   3 \tab Tukey \tab 1/3 \cr
#'   4 \tab Gumbel \tab 1 \cr
#'   5 \tab Hazen \tab 1/2 \cr
#'   6 \tab Cunnane \tab 2/5
#' }
#' 
#' @return
#' Matrix of probabilities per category
#' 
#' @references 
#' Wilks, D.S. (2011). Statistical methods in the atmospheric sciences (Third Edition). 
#' Academic press. 
#' 
#' @examples
#' tm <- toymodel()
#' 
#' ## convert to tercile forecasts (only display first forecast and obs)
#' count2prob(convert2prob(tm$fcst, prob=1:2/3))[1,]
#' count2prob(convert2prob(tm$obs, prob=1:2/3))[1,]
#' 
#' 
#' @seealso \code{\link{convert2prob}} for conversion of continuous forecasts to ensemble counts
#' 
#' @keywords utilities
#' @export
count2prob <- function(x, type=3){
  stopifnot(is.matrix(x))
  stopifnot(any(!is.na(x)))
  stopifnot(type %in% 1:6)
  stopifnot(rowSums(x) %% 1 == 0)

  if (all(rowSums(x) == 1)){
    xout <- x
  } else {
    ## select parameters to determine plotting position
    a <- c(0, 0.3, 1/3, 1, 1/2, 2/5)[type]
    n <- rowSums(x) + 1
    
    xout <- (x + 1 - a) / (n + 1 - 2*a)    
  }

  return(xout)
}