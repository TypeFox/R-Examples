################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015
##
##   This file is part of the R package reda.
##
##   The R package reda is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package reda is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################


#' Recurrent Event Data Analysis
#' 
#' The package \pkg{reda} mainly provides function \code{\link{rateReg}}
#' to fit parametric gamma frailty model with spline or piecewise constant
#' baseline rate function.
#' Another main function \code{\link{mcf}} computes and plots
#' the parametric mean cumulative function (MCF) from a fitted model
#' as well as the nonparametric sample MCF (Nelson-Aelson estimator)
#' for recurrent event data.
#'
#' See vignettes for introduction and demonstration.  
#'
#' @importFrom methods setClass setGeneric setMethod new validObject
#' @docType package
#' @name reda-package
NULL


