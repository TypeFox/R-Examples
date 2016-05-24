# EnsErrorss.R Ensemble Mean Error Skill Scores
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

#' @name EnsErrorss
#'   
#' @title Ensemble Mean Error Skill scores
#' 
#' @description Computes various ensemble mean error skill scores.
#' \code{EnsMess} computes the mean error, \code{EnsMaess} the mean absolute
#' error, \code{EnsMsess} the mean squared error, and \code{EnsRmsess} the
#' square root of the mean squared error (for consistency with the veri
#' package).
#' 
#' @param ens n x k matrix of n forecasts from k ensemble members
#' @param ens.ref n x l matrix of m reference forecasts from l ensemble members
#' @param obs n verifying observations
#' @param type specifying what error metric to compute, one of [me, mae, mse,
#'   rmse]
#'   
#' @examples
#' tm <- toymodel()
#' 
#' ## compute RMSE skill score against reference forecast with a bias of +2
#' EnsErrorss(ens=tm$fcst, ens.ref=tm$fcst + 2, obs=tm$obs, type='rmse')
#' 
#' ## compute skill score using veriApply
#' veriApply("EnsRmsess", fcst=tm$fcst, obs=tm$obs, fcst.ref=tm$fcst + 2)
#' 
#' 
#' @seealso \code{\link{veriApply}}, \code{\link{EnsError}}
#'   
#' @export
EnsErrorss <- function(ens, ens.ref, obs, type){
  stopifnot(is.matrix(ens), is.matrix(ens.ref), 
            is.vector(obs), length(obs) == nrow(ens),
            length(obs) == nrow(ens.ref))
  xmask <- apply(!is.na(ens), 1, any) & !is.na(obs) & apply(!is.na(ens.ref), 1, any)
  if (all(!xmask)) {
    xout <- NA
  } else {
    xout <- 1 - EnsError(ens, obs, type) / EnsError(ens.ref, obs, type)
  }
  return(xout)
}

#' @rdname EnsErrorss
#' @export
EnsMess <- function(ens, ens.ref, obs){
  EnsErrorss(ens=ens, ens.ref=ens.ref, obs=obs, type='me')
}

#' @rdname EnsErrorss
#' @export
EnsMaess <- function(ens, ens.ref, obs){
  EnsErrorss(ens=ens, ens.ref=ens.ref, obs=obs, type='mae')
}

#' @rdname EnsErrorss
#' @export
EnsMsess <- function(ens, ens.ref, obs){
  EnsErrorss(ens=ens, ens.ref=ens.ref, obs=obs, type='mse')
}

#' @rdname EnsErrorss
#' @export
EnsRmsess <- function(ens, ens.ref, obs){
  EnsErrorss(ens=ens, ens.ref=ens.ref, obs=obs, type='rmse')
}
