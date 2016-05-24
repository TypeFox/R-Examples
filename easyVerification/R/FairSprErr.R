# FairSprErr.R Spread to Error Ratio
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

#' Fair Spread to Error Ratio
#' 
#' Compute the spread to error ratio (\code{SPR}) for probabilistic forecasts -
#' not unlike the functions in SpecsVerification. \code{SPR > 1} indicates
#' overdispersion (underconfidence), whereas \code{SPR < 1} indicates
#' overconfidence in the forecasts.
#' 
#' @param ens n x k matrix of n forecasts for k ensemble members
#' @param obs vector with n verifying observations
#'   
#' @details Here we define the spread-error rate as the square root of the ratio
#'   of mean ensemble variance to the mean squared error of the ensemble mean 
#'   with the verifying observations. We inflate the intra ensemble sample 
#'   variance to account for the finite ensemble size as in Weigel (2011).
#'   
#'   
#' @references Weigel, A.P. (2012). Ensemble forecasts. Forecast Verification: A
#'   Practitioner's Guide in Atmospheric Science, Second Edition, 141-166.
#'   
#' @seealso \code{\link{veriApply}}, \code{\link{FairSprErr}}
#'   
#' @examples
#' tm <- toymodel()
#' FairSprErr(tm$fcst, tm$obs)
#' 
#' ## compute spread to error ratio using veriApply
#' veriApply('FairSprErr', fcst=tm$fcst, obs=tm$obs)
#' 
#' ## compare with 'unfair' spread to error ratio
#' veriApply("EnsSprErr", fcst=tm$fcst, obs=tm$obs)
#' 
#' @export
FairSprErr <- function(ens, obs){
  stopifnot(is.matrix(ens), is.vector(obs), nrow(ens) == length(obs))
  
  xmask <- apply(!is.na(ens), 1, any) & !is.na(obs)
  nens <- ncol(ens)
  spread <- mean(apply(ens[xmask,,drop=F], 1, sd, na.rm=T)**2, na.rm=T)
  error <- mean((obs - rowMeans(ens))**2, na.rm=T)
  
  return(sqrt((nens + 1) / nens * spread/error))
}