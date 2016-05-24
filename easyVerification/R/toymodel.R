# toymodel.R Create Example Forecast-observation Pairs
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

#' @name toymodel
#' @aliases toyarray
#' 
#' @title Create Example Forecast-Observation Pairs
#' 
#' @description
#' This toy model lets you create forecast-observation pairs with specified 
#' ensemble and forecast size, correlation skill, and overconfidence 
#' (underdispersion) for application with the verification functionality 
#' provided as part of the easyVerification package.
#' 
#' @param N number of forecast instances
#' @param nens number of ensemble members
#' @param alpha nominal correlation skill of forecasts
#' @param beta overconfidence parameter (see details)
#'   
#' @details The toy model is the TM2 model as introduced by Weigel and Bowler
#' (2009) with a slight modification to allow for forecasts with negative
#' correlation skill. In this toy model, the observations \eqn{x} and forecasts
#' \eqn{f_i} are defined as follows:
#' 
#' \deqn{x = \mu_x + \epsilon_x} 
#' \deqn{f_i = \alpha / |\alpha| \mu_x + \epsilon_{\beta} + \epsilon_i} 
#' where 
#' \deqn{\mu_x ~ N(0, \alpha^2)} 
#' \deqn{\epsilon_x ~ N(0, 1 - \alpha^2)} 
#' \deqn{\epsilon_{\beta} ~ N(0, \beta^2)} 
#' \deqn{\epsilon_i ~ N(0, 1 - \alpha^2 - \beta^2)} 
#' \deqn{\alpha^2 \le 1} 
#' \deqn{0 \le \beta \le 1 - \alpha^2}
#' 
#' @note This toy model is intended to provide example forecast observation
#' pairs and not to serve as a conceptual model to study real forecasts. For
#' models to do the latter, please refer to Siegert et al. (2015). 
#' 
#' @references A. Weigel and N. Bowler (2009). Comment on 'Can multi-model
#' combination really enhance the prediction skill of probabilistic ensemble
#' forecasts?'. \emph{Quarterly Journal of the Royal Meteorological Society},
#' 135, 535-539.
#' 
#' S. Siegert \emph{et al.} (2015). A Bayesian framework for verification and
#' recalibration of ensemble forecasts: How uncertain is NAO predictability?
#' Preprint on ArXiv, \url{http://arxiv.org/abs/1504.01933}.
#' 
#' 
#' @examples
#' ## compute the correlation for a toy forecast with default parameters
#' tm <- toyarray()
#' f.corr <- veriApply("EnsCorr", fcst=tm$fcst, obs=tm$obs)
#'  
#' @keywords utilities
#' @export
#' 
toymodel <- function(N=35, nens=51, alpha=0.5, beta=0){
  stopifnot(beta <= 1, beta >= 0, alpha**2 <= 1)
  if (alpha**2 + beta**2 > 1){
    stop("alpha^2 + beta^2 larger than 1")
  }
  signal <- rnorm(N, sd=sqrt(alpha**2))
  obs <- signal + rnorm(N, sd=sqrt(1 - alpha**2))
  fcst <- sign(alpha)*signal + rnorm(N, sd=sqrt(beta**2)) + 
    array(rnorm(N*nens, sd=sqrt(1 - alpha**2 - beta**2)), c(N, nens))
  return(list(obs=obs, fcst=fcst))
}

#' @rdname toymodel
#'
#' @param dims independent (e.g. spatial) dimensions for the toy model
#' @param ... additional arguments passed to \code{\link{toymodel}}
#' @export
toyarray <- function(dims=c(10,5), ...){
  tm <- sapply(rep(1, prod(dims)), function(x) toymodel(...))
  fcst <- sapply(tm['fcst', ], function(x) x, simplify='array')
  obs <- sapply(tm['obs', ], function(x) x)
  ## rearrange fcst and obs to put forecast and ensemble dimensions last
  fcst <- array(aperm(fcst, c(3,1,2)), c(dims, dim(fcst)[1:2]))
  obs <- array(t(obs), c(dims, nrow(obs)))
  return(list(fcst=fcst, obs=obs))
}