#
# file: plsr.mcSimulation.R
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
#' @include mcSimulation.R
NULL
##############################################################################################
# plsr.mcSimulation(object, ...)
##############################################################################################
#' Partial Least Squares Regression (PLSR) of Monte Carlo simulation results.
#' 
#' Perform a Partial Least Squares Regression (PLSR) of Monte Carlo simulation results.
#' @param object An object of class \code{mcSimulation}.
#' @param resultName \code{character}: indicating the name of the component of
#'   the simulation function (\code{model_function}) whose results histogram
#'   shall be generated. If \code{model_function} is single valued, no name
#'   needs to be supplied. Otherwise, one valid name has to be specified.
#'   Defaults to \code{NULL}.
#' @param variables.x \code{character} or \code{character} vector: Names of the 
#'   components of the input variables to the simulation function, i.e. the
#'   names of the variables in the input \code{estimate} which random sampling
#'   results shall be displayed. Defaults to all components.
#' @inheritParams pls::plsr
#' @param ... further arguments to be passed to \code{\link[pls:plsr]{plsr}}.
#' @return An object of class \code{\link[pls:mvr]{mvr}}.
#' @seealso \code{\link{mcSimulation}}, \code{\link[pls:plsr]{plsr}}, 
#' \code{\link[pls:summary.mvr]{summary.mvr}}, \code{\link[pls:biplot.mvr]{biplot.mvr}}, 
#' \code{\link[pls:coef.mvr]{coef.mvr}}, \code{\link[pls:plot.mvr]{plot.mvr}},
#' @export
plsr.mcSimulation <- function(object,
                              resultName=NULL,
                              variables.x=names(object$x),
                              method="oscorespls", scale=TRUE, ncomp=2,
                              ...
){
  # Check preconditions:
  ## Namespace requirements:
  if (!requireNamespace("pls", quietly = TRUE)) 
    stop("Package \"pls\" needed. Please install it.",
         call. = FALSE)
  # Prepare the dependent variable:
  if( is.list(object$y) ){
    if( !is.null(resultName) ){
      y<-object$y[[resultName]]
    } else {
      if(length(names(object$y))==1){
        y<-unlist(object$y)
      }
      else 
        stop("No component of the model function chosen!")
    } 
  } else { 
    y<-object$y
  }
  # Prepare the independent variables:
  x<-as.matrix((object$x)[variables.x])
  x<-x[,which(apply(x,2,sd)>0)]
  # PLSR analysis:
  plsrResults<-pls::plsr(y~x, method=method, scale=scale, ncomp=ncomp, ...) 
  
  # Return the PLSR results
  plsrResults
}
