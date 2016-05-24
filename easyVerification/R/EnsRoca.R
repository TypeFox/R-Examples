# EnsRoca.R Area Under the ROC Curve
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

#' Area Under the ROC Curve
#' 
#' Computes the area under the ROC curve given the observations
#' 
#' @param ens n x j matrix of n probability forecasts for j categories
#' @param obs n x j matrix of occurence of n verifying observations in j categories
#' 
#' @examples
#' tm <- toymodel()
#' 
#' ## compute ROC area for tercile forecasts using veriApply
#' veriApply("EnsRoca", fcst=tm$fcst, obs=tm$obs, prob=1:2/3)
#' 
#' @seealso \code{\link{veriApply}}, \code{\link{EnsRocss}}
#' 
#' @export
EnsRoca <- function(ens, obs){
  stopifnot(is.matrix(ens), is.matrix(obs), length(obs) == length(ens))
  roc.area <- EnsRocaCpp(ens, obs)
  roc.area <- as.list(roc.area)
  names(roc.area) <- paste0('cat', seq(along=roc.area))
  return(roc.area)
}

oldEnsRoca <- function(ens, obs){
  stopifnot(is.matrix(ens), is.matrix(obs), length(obs) == length(ens))
  rs <- rowSums(ens)
  if (any(rs != 1)){
    ## convert number of occurences to probabilities
    ens <- ens / rs
  }
  n.event <- colSums(obs)
  n.total <- nrow(obs)
  ens.rank <- apply(ens, 2, rank)
  mean.rank <- colSums(ens.rank * obs) / n.event
  roc.area <- (mean.rank - (n.event + 1)/2) / (n.total - n.event)
  roc.area[n.event == 0] <- NA
  roc.area <- as.list(roc.area)
  names(roc.area) <- paste0('cat', seq(along=roc.area))
  return(roc.area)
}