# Ens2AFC.R Generalized Discrimination Score
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

#' @name Ens2AFC
#'   
#' @title Generalized Discrimination Score
#' 
#' @description Computes the generalized discrimination score for ensemble
#' forecasts after (Weigel and Mason, 2011).
#' 
#' @param ens n x m matrix of n forecasts for m ensemble members
#' @param obs vector of n verifying observations
#' @param ... additional arguments not used in function (for compatibility)
#'   
#' @details This function computes the generalized discrimination score for
#' ensemble forecasts with continuous observations as described in Weigel and
#' Mason (2011).
#' 
#' @references Weigel, A.P., and S.J. Mason (2011). The Generalized
#' Discrimination Score for Ensemble Forecasts. Monthly Weather Review, 139(9),
#' 3069-3074. doi:10.1175/MWR-D-10-05069.1
#' 
#' @examples
#' tm <- toymodel()
#' Ens2AFC(tm$fcst, tm$obs)
#' 
#' @seealso \code{\link{veriApply}}
#'   
#' @export
Ens2AFC <- function(ens, obs, ...){
  return(0.5*(1 + cor(rankEnsCpp(ens), obs, method='kendall', use='p')))
}

#' @rdname Ens2AFC
rank.ensembles <- function (ens) {
  nens = dim(ens)[2]
  n = dim(ens)[1]
  ranks = rep(1, n)
  for (i in 2:n) for (j in 1:(i - 1)) {
    ens.tmp.event = ens[i, ]
    ens.tmp.nonev = ens[j, ]
    rank.1 = rank(c(ens.tmp.event, ens.tmp.nonev))[1:nens]
    p.afc = (sum(rank.1) - nens * (nens + 1)/2)/nens^2
    if (p.afc > 0.5) 
      ranks[i] = ranks[i] + 1
    if (p.afc < 0.5) 
      ranks[j] = ranks[j] + 1
    if (p.afc == 0.5) {
      ranks[i] = ranks[i] + 0.5
      ranks[j] = ranks[j] + 0.5
    }
  }
  return(ranks)
}

# rank.ens <- function(ens){
#   nens <- ncol(ens)
#   n <- nrow(ens)
#   ens.list <- apply(ens, 1, list)
#   U <- array(0.5, rep(n,2))
#   U[,] <- (apply(cbind(ens[rep(1:n, n),], ens[rep(1:n, each=n),]), 1, function(x) sum(rank(x)[1:nens])) - nens*(nens + 1)/2)/nens**2
#   ranks <- 0.5 + apply(sign(U - 0.5)/2 + 0.5, 2, sum)
#   return(ranks)
# }
