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


#' Simulated Sample Dataset for Demonstration
#'
#' A simulated data frame with covariates named 
#' \code{ID}, \code{time}, \code{event}, \code{group} and \code{x1},
#' where
#' \itemize{
#'     \item \code{ID}: Subjects identification;
#'     \item \code{time}: Event or censoring time;
#'     \item \code{event}: Event indicator, 1 = event, 0 = censored;
#'     \item \code{group}: Treatment group indicator;
#'     \item \code{x1}: Continuous variable. 
#' }
#'
#' @details
#' The sample dataset is originally simulated by the thinning
#' method developed by Lewis and Shedler (1979) and
#' further processed for a better demonstration purpose.
#' See Fu et al. (2014) for details also.
#' 
#' @docType data
#' @name simuDat
#' @format A data frame with 500 rows and 5 variables.
#' @references
#' Lewis, P. A., & Shedler, G. S. (1979). 
#' Simulation of nonhomogeneous Poisson processes by thinning.
#' \emph{Naval Research Logistics Quarterly}, 26(3), 403--413.
#' 
#' Fu, H., Luo, L., & Qu Y. (2014). Hypoglycemic Events Analysis via
#' Recurrent Time-to-Event (HEART) Models. 
#' \emph{Journal of biopharmaceutical statistics}, Epub 2014 Dec 1.
NULL

