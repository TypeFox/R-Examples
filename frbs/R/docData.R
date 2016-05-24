#############################################################################
#
#  This file is a part of the R package "frbs".
#
#  Author: Lala Septem Riza
#  Co-author: Christoph Bergmeir
#  Supervisors: Francisco Herrera Triguero and Jose Manuel Benitez
#  Copyright (c) DiCITS Lab, Sci2s group, DECSAI, University of Granada.
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
#' The package includes embedded versions of the Mackey-Glass chaotic time series and the Gas Furnance dataset.
#'
#' \bold{Mackey-Glass chaotic time series}
#' 
#' The Mackey-Glass chaotic time series is defined by the following delayed differential equation:
#' 
#' \eqn{d_x(t) / d_t = (a * x(t - \tau) / (1 + x(t - \tau) ^ 10)) - b * x(t)}
#' 
#' For this dataset, we generated 1000 samples, with input parameters as follows:
#' \itemize{
#' \item \eqn{a = 0.2}
#' \item \eqn{b = 0.1}
#' \item \eqn{\tau = 17}
#' \item \eqn{x_0 = 1.2}
#' \item \eqn{d_t = 1}
#' }
#' 
#' The dataset is embedded in the following way: 
#'
#' input variables: \eqn{x(t - 18)}, \eqn{x(t - 12)}, \eqn{x(t - 6)}, \eqn{x(t)}
#'
#' output variable: \eqn{x(t + 6)}
#'
#' \bold{Gas Furnance dataset}
#' 
#' The Gas Furnance dataset is taken from Box and Jenkins. It consists of 292 consecutive 
#' values of methane at time \eqn{(t - 4)}, and the CO2 produced in a furnance at time \eqn{(t - 1)} as input 
#' variables, with the produced CO2 at time \eqn{(t)} as an output variable. So, each training data 
#' point consists of \eqn{[u(t - 4), y(t - 1), y(t)]}, where \eqn{u} is methane and \eqn{y} is CO2.
#'
#' @title Data set of the package
#' @name frbsData
#' @docType data
#' @references 
#' G. E. P. Box and G. M. Jenkins, "Time series analysis, forecasting and control", San Fransisco, CA: Holden Day (1970).
#' 
#' M. Mackey and L. Glass, "Oscillation and chaos in physiological control systems", Science, vol. 197, pp. 287 - 289 (1977).
#' 
#' @keywords data
NULL