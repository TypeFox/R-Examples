#
# file: rtnorm90ci.R
#
# This file is part of the R-package decisionSupport
# 
# Authors: 
#   Lutz Göhring <lutz.goehring@gmx.de>
#   Eike Luedeling (ICRAF) <eike@eikeluedeling.com>
#
# Copyright (C) 2015 World Agroforestry Centre (ICRAF) 
#	http://www.worldagroforestry.org
# 
# The R-package decisionSupport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# The R-package decisionSupport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with the R-package decisionSupport.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################################
#' @include paramtnormci_numeric.R
#' @include paramtnormci_fit.R
NULL
##############################################################################################
# rtnorm90ci(n, lower, upper, relativeTolerance)
##############################################################################################
#' 90\%-confidence interval based truncated normal random number generation.
#' 
#' \code{rtnorm90ci} generates truncated normal random numbers based on the 90\% confidence interval
#' calculating the distribution parameter numerically from the  90\%-confidence interval or via a 
#' fit on the 90\%-confidence interval. The fit might include the median or not.
#' @param n Number of generated observations.
#' @param ci \code{numeric} 2-dimensional vector; lower, i.e \code{ci[[1]]}, and upper bound, i.e
#'   \code{ci[[2]]}, of the  90\%-confidence interval.
#' @param median if \code{NULL}: truncated normal is fitted only to lower and upper value of the 
#'   confidence interval; if \code{numeric}: truncated normal is fitted on the confidence interval 
#'   and the median simultaneously. For details cf. below. This option is only relevant if 
#'   \code{method="fit"}.
#' @param lowerTrunc \code{numeric}; lower truncation point of the distribution (>= \code{-Inf}).
#' @param upperTrunc \code{numeric}; upper truncation point of the distribution (<= \code{Inf}).
#' @param method method used to determine the parameters of the truncated normal; possible methods 
#'   are \code{"numeric"} (the default) and \code{"fit"}.
#' @param relativeTolerance \code{numeric}; the relative tolerance level of deviation of the
#'   generated confidence interval from the specified interval. If this deviation is greater than
#'   \code{relativeTolerance} a warning is given.
#' @param ... further parameters to be passed to \code{\link{paramtnormci_numeric}} or 
#'   \code{\link{paramtnormci_fit}}, respectively.
#' @details
#' \code{method="numeric"} is implemented by \code{\link{paramtnormci_numeric}} and 
#' \code{method="fit"} by \code{\link{paramtnormci_fit}}.
#' @seealso For the implementation of \code{method="numeric"}: \code{\link{paramtnormci_numeric}}; 
#'  for the implementation of \code{method="fit"}: \code{\link{paramtnormci_fit}}.
#' @export
rtnorm90ci <- function(n, ci, median=mean(ci), lowerTrunc=-Inf, upperTrunc=Inf, method="numeric", relativeTolerance=0.05,...){
  # Constants:
  p=c(0.05, 0.95)
  # Create output vector for the random numbers to be generated
  x<-vector(length=n)
  # Calculate mean and sd corresponding to confidence interval according to method:
  if( method=="numeric" )
    param<-paramtnormci_numeric(p=p, ci=ci, lowerTrunc=lowerTrunc, upperTrunc=upperTrunc, 
                                relativeTolerance=relativeTolerance)
  else if(method=="fit")
    param<-paramtnormci_fit(p=p, ci=ci, median=median, lowerTrunc=lowerTrunc, upperTrunc=upperTrunc, 
                            relativeTolerance=relativeTolerance)
  else
    stop("method=\"", method, "\" does not exist!")
  # Generate the random numbers:
  x<-msm::rtnorm(n=n,
                 mean=param$mean,
                 sd=param$sd,
                 lower=lowerTrunc,
                 upper=upperTrunc)
  #Return the generated positive normal random numbers:
  x
}
##############################################################################################
# rposnorm90ci(n, lower, upper, relativeTolerance)
##############################################################################################
#' 90\%-confidence interval based positive normal random number generation.
#' 
#' \code{rposnorm90ci} generates positive normal random numbers based on the 90\% confidence interval. 
#' It is a wrapper function for \code{rtnorm90ci}.
#' @rdname rtnorm90ci
#' @param lower \code{numeric}; lower bound of the 90\% confidence interval.
#' @param upper \code{numeric}; upper bound of the 90\% confidence interval. 
#' @details 
#' Positive normal random number generation: a positive normal distribution
#' is a truncated normal distribution with lower truncation point equal to zero and upper truncation
#' is infinity. \code{rposnorm90ci} implements this as a wrapper function for \ifelse{latex}{\cr}{ }
#' \code{rtnorm90ci(n, c(lower,upper), median, lowerTrunc=0, upperTrunc=Inf, method, relativeTolerance,...)}.
#' @export
rposnorm90ci <- function(n, lower, median=mean(c(lower,upper)), upper, method="numeric", relativeTolerance=0.05,...){
 lowerTrunc<-0
 upperTrunc<-Inf
  # Generate the random numbers:
  x<-rtnorm90ci(n=n,
                ci=(c(lower, upper)),
                median=median,
                lowerTrunc=lowerTrunc,
                upperTrunc=upperTrunc,
                method=method, 
                relativeTolerance=relativeTolerance,
                ...)
  #Return the generated positive normal random numbers:
  x
}
##############################################################################################
# rtnorm_0_1_90ci(n, lower, upper, relativeTolerance)
##############################################################################################
#' 90\%-confidence interval based normal random number generation truncated to \eqn{[0,1]}.
#' 
#' \code{rtnorm_0_1_90ci} generates normal random numbers truncated to \eqn{[0,1]} based on the 
#' 90\% confidence interval. It is a wrapper function for \code{rtnorm90ci}.
#' @rdname rtnorm90ci
#' @details 
#' 0-1-(truncated) normal random number generation: a 0-1-normal distribution
#' is a truncated normal distribution with lower truncation point equal to zero and upper truncation
#' equal to 1. \code{rtnorm_0_1_90ci} implements this as a wrapper function for \ifelse{latex}{\cr}{ }
#' \code{rtnorm90ci(n, c(lower,upper), median, lowerTrunc=0, upperTrunc=1, method, relativeTolerance,...)}.
#' @export
rtnorm_0_1_90ci <- function(n, lower, median=mean(c(lower,upper)), upper, method="numeric", relativeTolerance=0.05,...){
  lowerTrunc<-0
  upperTrunc<-1
  # Generate the random numbers:
  x<-rtnorm90ci(n=n,
                ci=(c(lower, upper)),
                median=median,
                lowerTrunc=lowerTrunc,
                upperTrunc=upperTrunc,
                method=method, 
                relativeTolerance=relativeTolerance,
                ...)
  #Return the generated positive normal random numbers:
  x
}

