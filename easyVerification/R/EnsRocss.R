# EnsRocss.R ROC Area Skill Score
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

#' Skill Score for Area under the ROC Curve
#' 
#' Computes the skill score for the area under the ROC curve compared to an 
#' arbitrary reference forecast (generally climatological forecast).
#' 
#' @param ens n x j matrix of n probability forecasts for j categories
#' @param ens.ref n x j matrix of reference forecast for j categories
#' @param obs n x j matrix of occurence of n verifying observations in j 
#'   categories
#'   
#' @return a list with the ROC area skill score for forecast category \code{i} 
#'   in \code{cati} and the standard deviation of this skill score for category 
#'   \code{i} in \code{cati.sigma} if a reference forecast with zero association
#'   is used (see details).
#'   
#' @details For the traditional ROC area skill score where the reference 
#'   forecast has zero association with the observations, the standard error 
#'   \eqn{\sigma} of the ROC area skill score is given by the following formula 
#'   after Broecker (2012).
#'   
#'   \deqn{\sigma^2 = \frac{1}{3} \left(\frac{1}{N_0} + \frac{1}{N_1} + 
#'   \frac{1}{N_0 N_1} \right)}{\sigma^2 = 1/3 (1/N0 + 1/N1 + 1/(N0 N1))}
#'   
#'   Where \eqn{\sigma} is the standard error, \eqn{N_1}{N1} the number of 
#'   events, and \eqn{N_0}{N0} the number of non-events in category \code{i}. 
#'   Please note the factor 2 difference to the formulation of the standard 
#'   error for the ROC area in the original manuscript due to the conversion of
#'   the ROC area to the ROC area skill score.
#'   
#' @references Br\"ocker, J. (2012). Probability forecasts. Forecast 
#'   Verification: A Practitioner's Guide in Atmospheric Science, Second 
#'   Edition, 119-139.
#'   
#' @examples
#' tm <- toymodel()
#' 
#' ## compute ROC skill score for forecasts of x <= 0, 0 <= x < 1, and x > 1
#' ## skill score is computed using climatological forecast as reference
#' veriApply("EnsRocss", tm$fcst, tm$obs, threshold=c(0,1))
#'   
#' @seealso \code{\link{veriApply}}, \code{\link{EnsRoca}}
#'   
#' @export
EnsRocss <- function(ens, ens.ref, obs){
  roc.area <- EnsRoca(ens, obs)
  if (all(ens.ref == rep(ens.ref[1,], each=nrow(ens.ref)))){
    roc.out <- lapply(roc.area, function(x) 2*x - 1)   
    ## compute sigma
    N1 <- apply(obs, 2, sum)
    N0 <- nrow(obs) - N1
    roc.sigma <- 2 * sqrt(1/12 * (1/N0 + 1/N1 + 1/(N0*N1)))
    roc.sigma[N1 == 0] <- NA
    roc.sigma <- as.list(roc.sigma)
    names(roc.sigma) <- paste0('cat', seq(along=roc.sigma), '.sigma')
    roc.out <- c(roc.out, roc.sigma)
  } else {
    roc.ref <- EnsRoca(ens.ref, obs)
    roc.out <- Map(function(x,y) x/y - 1, x=roc.area, y=roc.ref)
  }
  
  return(roc.out)
}