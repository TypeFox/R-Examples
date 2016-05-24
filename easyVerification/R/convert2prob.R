# convert2prob.R Convert to Category Forecast
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

#' Convert to Probability / Category Forecast
#' 
#' Converts the continuous ensemble forecast to counts of ensemble members per
#' category. The categories can be defined relative to the ensemble distribution
#' (using \code{prob}) or relative to absolute values for the category
#' thresholds (using \code{threshold}, see details).
#' 
#' @param x input vector or matrix
#' @param prob thresholds for categorical forecasts (defaults to NULL)
#' @param threshold absolute thresholds for categorical forecasts (defaults to
#'   NULL)
#' @param multi.model logical, are we dealing with initial condition (the
#'   default) or multi-model ensembles (see details)?
#'   
#' @details In case both \code{prob} and \code{threshold} are set to
#' \code{NULL}, the function returns the input \code{x} without modification. If
#' \code{prob} is set, a matrix with the number of occurences per class for a
#' given quantile of the full distribution (e.g. temperature above/below the
#' median). If \code{threshold} is set, the classes are defined based on the
#' absolute value (e.g. temperature above/below 13 deg. C). Multiple classes are
#' supported.
#' 
#' If \code{multi.model = TRUE}, the relative thresholds supplied by \code{prob}
#' are ensemble member specific, i.e. are estimated for each ensemble member 
#' separately. This is in particular applicable for multi-model ensembles with 
#' model dependent biases.
#' 
#' @return Matrix of occurences per class (i.e. the number of ensemble members
#' per class, or an indicator for the observations)
#' 
#' @examples
#' tm <- toymodel()
#' 
#' ## convert to tercile forecasts (only display first forecast and obs)
#' convert2prob(tm$fcst, prob=1:2/3)[1,]
#' convert2prob(tm$obs, prob=1:2/3)[1,]
#' 
#' ## convert to category forecasts (smaller and larger than 1)
#' convert2prob(tm$fcst, threshold=1)[1,]
#' convert2prob(tm$obs, threshold=1)[1,]
#' 
#' @seealso \code{\link{veriApply}}
#'   
#' @keywords utilities
#' @export
convert2prob <- function(x, prob=NULL, threshold=NULL, multi.model=FALSE){
  stopifnot(is.vector(x) | is.matrix(x))
  stopifnot(any(!is.na(x)))
  if (!is.null(prob) & !is.null(threshold)){
    stop('Both probability and absolute thresholds provided')
  } 
  ## convert probability to absolute threshold
  if (is.numeric(prob)){
    if (multi.model){
      threshold <- apply(x, 2, quantile, sort(prob), na.rm=T, type=8)
    } else {
      threshold <- quantile(x, sort(prob), na.rm=T, type=8)      
    }
  }
  ## compute occurence per class
  if (is.numeric(threshold)){
    nens <- ncol(as.matrix(x))
    nclass <- nrow(as.matrix(threshold)) + 1
    if (multi.model & !is.null(prob) & nens > 1){
      xtmp <- array(sapply(1:nens, function(i) 
        apply(sapply(threshold[,i], function(y) x[,i] > y), 1, sum)), 
        dim(as.matrix(x))) + 1
    } else {
      xtmp <- array(apply(sapply(sort(threshold), function(y) c(x) > y), 1, sum), dim(as.matrix(x))) + 1
    }
    xout <- t(apply(xtmp, 1, tabulate, nbins=nclass))      
    xout[apply(as.matrix(is.na(x)), 1, any),] <- NA
  } else {
    xout <- x
  }
  return(xout)
}