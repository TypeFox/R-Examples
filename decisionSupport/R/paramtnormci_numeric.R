#
# file: paramtnormci_numeric.R
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
# paramtnormci_numeric(p, ci, lowerTrunc, upperTrunc, relativeTolerance)
##############################################################################################
#' Return parameters of truncated normal distribution based on a confidence interval.
#'
#' This function calculates the distribution parameters, i.e. \code{mean} and \code{sd}, of a
#' truncated normal distribution from an arbitrary confidence interval.
#' @param p \code{numeric} 2-dimensional vector; probabilities of lower and upper bound of the
#'   corresponding confidence interval.
#' @param ci \code{numeric} 2-dimensional vector; lower, i.e \code{ci[[1]]}, and upper bound, i.e
#'   \code{ci[[2]]}, of the  confidence interval.
#' @param lowerTrunc \code{numeric}; lower truncation point of the distribution (>= \code{-Inf}).
#' @param upperTrunc \code{numeric}; upper truncation point of the distribution (<= \code{Inf}).
#' @param relativeTolerance \code{numeric}; the relative tolerance level of deviation of the
#'   generated confidence interval from the specified interval. If this deviation is greater than
#'   \code{relativeTolerance} a warning is given.
#' @param rootMethod \code{character}; if \code{="probability"} the equation defining the parameters \code{mean} and
#'   \code{sd} is the difference between calculated and given probabilities of the confidence
#'   interval; if \code{="quantile"} the equation defining the parameters is the difference between
#'   calculated and given upper and lower value of the confidence interval.
#' @param ... Further parameters passed to \code{\link[nleqslv]{nleqslv}}.
#' @return A list with elements \code{mean} and \code{sd}, i.e. the parameters of the underlying
#'   normal distribution.
#' @details For details of the truncated normal distribution see \code{\link[msm]{tnorm}}.
#'
#' @seealso \code{\link[msm]{tnorm}}, \code{\link[nleqslv]{nleqslv}}
#' @export
paramtnormci_numeric <- function(p, ci, lowerTrunc=-Inf, upperTrunc=Inf, relativeTolerance=0.05,
                                 rootMethod="probability", ...){
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
  
  # Initialize the root finding:
  mean_init <- mean(ci)
  sd_init<- (mean_init - ci[["lower"]])/c_0.95
  
  if ( rootMethod=="quantile"){
    # Function defined by the difference between the target confidence values and the calculated
    # confidence values for certain values of the parameters mean and sd. Thus this function defines
    # mean and sd by f_calc(x) = 0, (x[1]:=mean, x[2]:=sd):
    f_calc <-function(x){
      msm::qtnorm(p=p, mean=x[1], sd=x[2], lower=lowerTrunc, upper=upperTrunc) - ci
    }
    # Fall back function for f_calc by random sampling simulation, i.e. function defined by the
    # difference between the target confidence values and the simulated confidence values for
    # certain values of the parameters mean and sd. Thus this function defines mean and sd by
    # f_calc(x) = 0, (x[1]:=mean, x[2]:=sd):
    f_sim<-function(x){
      n<-100*as.integer(1/(relativeTolerance*relativeTolerance))
      r<- msm::rtnorm(n=n, mean=x[1], sd=x[2], lower=lowerTrunc, upper=upperTrunc)
      
      quantile(x=r,probs=p) - ci
    }
    
  } else if( rootMethod=="probability"){
    # Function defined by the difference between confidence probabilities p and the calculated
    # probability for certain values of the parameters mean and sd. Thus this function defines
    # mean and sd by f_calc(x) = 0, (x[1]:=mean, x[2]:=sd):
    f_calc <-function(x){
      y <- msm::ptnorm(q=ci, mean=x[1], sd=x[2], lower=lowerTrunc, upper=upperTrunc) - p
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
      
      length(r[ r<= ci ])/n - p
    }
    
  } else
    stop("No root finding method chosen.")
  # Function wrapping f_calc and f_sim and thus defining mean and sd by f(x) = 0
  # (x[1]:=mean, x[2]:=sd):
  f <- function(x){
    tryCatch(f_calc(x=x),
             error=function(e) f_sim(x=x)
    )
  }
  
  # The root of f are mean and sd:
  #	x_0<-nleqslv::nleqslv(x=c(mean_init, sd_init), fn=f, control=list(maxit=10000))
  x_0<-nleqslv::nleqslv(x=c(mean_init, sd_init), fn=f, ...)
  mean<-x_0$x[1]
  sd<-x_0$x[2]
  
  
  # Check postcondition:
  tryCatch( ci_calc<- msm::qtnorm(p=p, mean=mean, sd=sd, lower=lowerTrunc, upper=upperTrunc),
            error=function(e){
              n<-100*as.integer(1/(relativeTolerance*relativeTolerance))
              r<- msm::rtnorm(n=n, mean=mean, sd=sd, lower=lowerTrunc, upper=upperTrunc)
              ci_calc<- quantile(x=r,probs=p)
            }
  )
  p_calc<-msm::ptnorm(q=ci, mean=mean, sd=sd, lower=lowerTrunc, upper=upperTrunc)
  for( j in seq(along=p) ){
    scale <- if( p[[j]] > 0 ) p[[j]] else NULL
    if( !isTRUE( msg<-all.equal(p[[j]], p_calc[[j]],  scale=scale, tolerance=relativeTolerance) ) ){
        warning("Calculated value of ", 100*p[[j]], "%-quantile: ", ci_calc[[j]], "\n  ",
                "Target value of ", 100*p[[j]], "%-quantile:     ", ci[[j]],   "\n  ",
                "Calculated cumulative probability at value ", ci[[j]], " : ", p_calc[[j]], "\n  ",
                "Target  cumulative probability at value ", ci[[j]], " : ", p[[j]], "\n  ",
                msg)
    }
  }
  
  #Return the calculated parameters:
  list(mean=mean, sd=sd)
}

