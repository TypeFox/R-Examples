# EnsError.R Ensemble Mean Error
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

#' @name EnsError
#'   
#' @title Ensemble Mean Error
#' 
#' @description Computes various ensemble mean error scores. \code{EnsMe}
#' computes the mean error, \code{EnsMae} the mean absolute error, \code{EnsMse}
#' the mean squared error, and \code{EnsRmse} the square root of the mean
#' squared error (for consistency with the veri package).
#' 
#' @param ens n x k matrix of n forecasts from k ensemble members
#' @param obs n verifying observations
#' @param type specifying what error metric to compute, one of [me, mae, mse,
#'   rmse]
#'   
#' @examples
#' #forecast and observations
#' tm <- toymodel()
#' 
#' # compute the mean bias
#' EnsError(tm$fcst, tm$obs, type='me')
#' # equivalently
#' EnsMe(tm$fcst, tm$obs)
#' 
#' @seealso \code{\link{veriApply}}, \code{\link{EnsErrorss}}
#'   
#' @export
EnsError <- function(ens, obs, type){
  stopifnot(is.matrix(ens), is.vector(obs), length(obs) == nrow(ens))
  xmask <- apply(!is.na(ens), 1, any) & !is.na(obs)
  if (all(!xmask)) {
    xout <- NA
  } else {
    error <- rowMeans(ens) - obs    
    if (type == 'me'){
      xout <- mean(error, na.rm=T)
    } else if (type == 'mae') {
      xout <- mean(abs(error), na.rm=T)
    } else if (type == 'mse') {
      xout <- mean(error**2, na.rm=T)
    } else if (type == 'rmse') {
      xout <- sqrt(mean(error**2, na.rm=T))
    }
  }
  return(xout)
}

#' @rdname EnsError
#' @export
EnsMe <- function(ens, obs){
  EnsError(ens=ens, obs=obs, type='me')
}

#' @rdname EnsError
#' @export
EnsMae <- function(ens, obs){
  EnsError(ens=ens, obs=obs, type='mae')
}

#' @rdname EnsError
#' @export
EnsMse <- function(ens, obs){
  EnsError(ens=ens, obs=obs, type='mse')
}

#' @rdname EnsError
#' @export
EnsRmse <- function(ens, obs){
  EnsError(ens=ens, obs=obs, type='rmse')
}
