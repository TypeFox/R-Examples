#
# file: paramtnormci_fit.R
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
##############################################################################################
# paramtnormci_fit(p, ci, median, lowerTrunc, upperTrunc, relativeTolerance)
##############################################################################################
#' Fit parameters of truncated normal distribution based on a confidence interval.
#'
#' This function fits the distribution parameters, i.e. \code{mean} and \code{sd}, of a truncated
#' normal distribution from an arbitrary confidence interval and, optionally, the median.
#' @param p \code{numeric} 2-dimensional vector; probabilities of upper and lower bound of the
#'   corresponding confidence interval.
#' @param ci \code{numeric} 2-dimensional vector; lower, i.e \code{ci[[1]]}, and upper bound, i.e
#'   \code{ci[[2]]}, of the  confidence interval.
#' @param  median if \code{NULL}: truncated normal is fitted only to lower and upper value of the
#'   confidence interval; if \code{numeric}: truncated normal is fitted on the confidence interval
#'   and the median simultaneously. For details cf. below.
#' @param lowerTrunc \code{numeric}; lower truncation point of the distribution (>= \code{-Inf}).
#' @param upperTrunc \code{numeric}; upper truncation point of the distribution (<= \code{Inf}).
#' @param relativeTolerance \code{numeric}; the relative tolerance level of deviation of the
#'   generated probability levels from the specified confidence interval. If the relative deviation
#'   is greater than \code{relativeTolerance} a warning is given.
#' @param fitMethod optimization method used in \code{\link[stats]{constrOptim}}.
#' @param ... further parameters to be passed to \code{\link[stats]{constrOptim}}.
#' @return A list with elements \code{mean} and \code{sd}, i.e. the parameters of the underlying
#'   normal distribution.
#' @details For details of the truncated normal distribution see \code{\link[msm]{tnorm}}.
#'
#'   The cumulative distribution of a truncated normal \eqn{F_{\mu, \sigma}}(x) gives the
#'   probability that a sampled value is less than \eqn{x}. This is equivalent to saying that for
#'   the vector of quantiles \ifelse{latex}{\eqn{q=(q_{p_1}, \ldots, q_{p_k})}}{\eqn{q=(q(p_1),
#'   \ldots, q(p_k))}} at the corresponding probabilities \eqn{p=(p_1, \ldots, p_k)} it holds that
#'     \deqn{p_i = F_{\mu, \sigma}(q_{p_i}),~i = 1, \ldots, k}{p_i = F_{\mu, \sigma}(q(p_i)), i = 1, \ldots k.}
#'   In the case of arbitrary postulated quantiles this system of equations might not have a
#'   solution in \eqn{\mu} and \eqn{\sigma}. A least squares fit leads to an approximate solution:
#'     \deqn{\sum_{i=1}^k (p_i - F_{\mu, \sigma}(q_{p_i}))^2 = \min}{ \sum_{i=1}^k (p_i - F_{\mu, \sigma}(q(p_i)))^2 = min}
#'   defines the parameters \eqn{\mu} and \eqn{\sigma} of the underlying normal distribution. This
#'   method solves this minimization problem for two cases:
#'    \enumerate{
#'      \item \code{ci[[1]] < median < ci[[2]]}: The parameters are fitted on the lower and upper value
#'        of the confidence interval and the median, formally:\cr
#'        \eqn{k=3}\cr
#'        \eqn{p_1}=\code{p[[1]]}, \eqn{p_2}=\code{0.5} and \eqn{p_3}=\code{p[[2]]};\cr
#'        \ifelse{latex}{\eqn{q_{p_1}}}{\eqn{q(p_1)}}=\code{ci[[1]]},
#'        \ifelse{latex}{\eqn{q_{0.5}}}{\eqn{q(0.5)}}=\code{median} and
#'        \ifelse{latex}{\eqn{q_{p_3}}}{\eqn{q(p_3)}}=\code{ci[[2]]}
#'      \item \code{median=NULL}: The parameters are fitted on the lower and upper value of the
#'        confidence interval only, formally:\cr
#'        \eqn{k=2}\cr
#'        \eqn{p_1}=\code{p[[1]]}, \eqn{p_2}=\code{p[[2]]};\cr
#'        \ifelse{latex}{\eqn{q_{p_1}}}{\eqn{q(p_1)}}=\code{ci[[1]]},
#'        \ifelse{latex}{\eqn{q_{p_2}}}{\eqn{q(p_2)}}=\code{ci[[2]]}
#'    }
#'    The \code{(p[[2]]-p[[1]])} - confidence interval must be symmetric in the sense that
#'    \code{p[[1]] + p[[2]] = 1}.
#'
#' @seealso \code{\link[msm]{tnorm}}, \code{\link[stats]{constrOptim}}
#' @export
paramtnormci_fit <- function(p, ci, median=mean(ci), lowerTrunc=-Inf, upperTrunc=Inf, relativeTolerance=0.05,
                             fitMethod="Nelder-Mead", ...){
  # Namespace requirements:
  requiredPackage<-"msm"
  if( !requireNamespace(requiredPackage, quietly = TRUE) )
    stop("Package \"",requiredPackage,"\" needed for truncated normal distributions. Please install it.",
         call. = FALSE)
  # Constants:
  # 95%-critical value of standard normal distribution (c_0.95=1.645):
  c_0.95=qnorm(0.95)
  # Check preconditions
  if ( is.null(p) || !all(!is.na(p)))
    stop("p must be supplied.")
  if ( is.null(ci) || !all(!is.na(ci)))
    stop("ci must be supplied.")
  if ( is.null(lowerTrunc) || is.null(upperTrunc) || is.na(lowerTrunc) || is.na(upperTrunc) )
    stop("lower and upper truncation points must be supplied.")
  if (length(p)!=2)
    stop("p must be of length 2.")
  if (length(ci)!=2)
    stop("ci must be of length 2.")
  # Prepare input variable: types
  p<-as.numeric(p)
  ci<-as.numeric(ci)
  lowerTrunc<-as.numeric(lowerTrunc)
  upperTrunc<-as.numeric(upperTrunc)
  if(p[[1]] >= p[[2]])
    stop("p[[1]] >= p[[2]]")
  if(ci[[1]] >= ci[[2]])
    stop("ci[[1]] >= ci[[2]]")
  names(p)<-c("lower", "upper")
  names(ci)<-c("lower", "upper")
  if ( !((lowerTrunc < ci[["lower"]] &&  ci[["upper"]] < upperTrunc))  )
    stop("ci is not a subset of [lowerTrunc, upperTrunc]!")
  if( sum(p) != 1)
    stop( "sum(p) != 1 " )
  # Quantile definition:
  if(!is.null(median) && !is.na(as.numeric(median))){
    median<-as.numeric(median)
    if ( !(( ci[["lower"]] < median &&  median < ci[["upper"]] ))  )
      stop("ci does not contain the median!")
    q<-c(ci[["lower"]], median, ci[["upper"]])
    p<-c(p[["lower"]], 0.5, p[["upper"]])
  } else {
    warning("median is not supplied; fitting the parameters of a truncated normal distribution only
            on the confidence interval might not lead to the desired distributions shape.")
    q<-c(ci[["lower"]], ci[["upper"]])
    p<-c(p[["lower"]],  p[["upper"]])
  }
  # Initialize the fit:
  mean_init <- if( !is.null(median) ) median else mean(ci)
  sd_init <- (mean(ci) - ci[["lower"]])/c_0.95

  # Function defined by the difference between confidence probabilities p and the calculated
  # probability for certain values of the parameters mean and sd. Thus this function defines
  # mean and sd by f_calc(x) = 0, (x[1]:=mean, x[2]:=sd):
  f_calc <-function(x){
    y <- msm::ptnorm(q=q, mean=x[1], sd=x[2], lower=lowerTrunc, upper=upperTrunc) - p
    # Produce error in case on NAs such that the function can be caught
    if (any(is.na(y))) stop ("NAs produced")
    y
  }
  # Fall back function for f_calc by random sampling simulation, i.e. function defined by the
  # difference between confidence probabilities p and the simulated probability for certain values
  # of the parameters mean and sd. Thus this function defines
  # mean and sd by f_calc(x) = 0, (x[1]:=mean, x[2]:=sd):
  f_sim<-function(x){
    n<-100*as.integer(1/(relativeTolerance*relativeTolerance))
    r<- msm::rtnorm(n=n, mean=x[1], sd=x[2], lower=lowerTrunc, upper=upperTrunc)
    y <- length(r[ r<= q ])/n - p
    y
  }
  # Function wrapping f_calc and f_sim and thus defining mean and sd by f(x) = 0
  # (x[1]:=mean, x[2]:=sd):
  f <- function(x){
    tryCatch(f_calc(x=x),
             error=function(e) f_sim(x=x)
    )
  }
  # The squared euclidean norm of the function f:
  g <- function(x) sum(f(x)*f(x))

  # The optimization constraint:
  # sd > 0 (for numerical stability)
  # is implemented as follows:
  # Constraints A x + B > 0
  # A_1 = ( 0  1 )
  # x_1 = mean
  # x_2 = sd
  # B_1 = 0
  A <- rbind(c( 0, 1))
  B <- rbind(c(  0 ))
  constraints<-list(ineqA=A, ineqB=B)

  # Fit the paramteters mean and sd:
  optimizationResult<-stats::constrOptim(theta=c(mean_init, sd_init),
                                         f=g,
                                         grad=NULL,
                                         ui=A,
                                         ci=-B,
                                         method=fitMethod, ...)
  # Save the fitted target parameters:
  mean<-optimizationResult$par[1]
  sd<-optimizationResult$par[2]
  # Check postcondition:
  tryCatch(q_calc <- msm::qtnorm(p, mean=mean, sd=sd, lower=lowerTrunc, upper=upperTrunc),
           error=function(e){
             n<-100*as.integer(1/(relativeTolerance*relativeTolerance))
             r<- msm::rtnorm(n=n, mean=mean, sd=sd, lower=lowerTrunc, upper=upperTrunc)
             q_calc <- quantile(x=r,probs=p)
           }
  )
  p_calc<-msm::ptnorm(q=q, mean=mean, sd=sd, lower=lowerTrunc, upper=upperTrunc)
  for( j in seq(along=p) ){
    scale <- if( p[[j]] > 0 ) p[[j]] else NULL
    if( !isTRUE( msg<-all.equal(p[[j]], p_calc[[j]],  scale=scale, tolerance=relativeTolerance) ) ){
      warning("Fitted value of ", 100*p[[j]], "%-quantile: ", q_calc[[j]], "\n  ",
              "Target value of ", 100*p[[j]], "%-quantile: ", q[[j]],   "\n  ",
              "Fitted cumulative probability at value ", q[[j]], " : ", p_calc[[j]], "\n  ",
              "Target cumulative probability at value ", q[[j]], " : ", p[[j]], "\n  ",
              msg)
    }
  }
  # Return the fitted parameters:
  list(mean=mean, sd=sd)
}


