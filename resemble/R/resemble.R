#' This is the version \code{1.2} (\code{'alma-de-coco'}) of the package. It 
#' implements a number of \code{R} functions useful for 
#' modeling complex visible and infrared spectra (\acronym{vis-IR}). 
#' The package includes functions for projecting spectral data 
#' onto orthogonal spaces, computing spectral dissimilarity
#' matrices, removing irrelevant spectra from a reference set, 
#' and modeling spectral data using memory-based learning.
#' 
#' The functions available for projecting the spectra are:
#' \itemize{
#'   \item{\code{\link{orthoProjection}}} 
#'   \item{\code{\link{pcProjection}}} 
#'   \item{\code{\link{plsProjection}}} 
#'   \item{\code{\link{predict.orthoProjection}}} 
#'   }
#' The functions available for computing similarity/dissimilarity 
#' matrices are:
#' \itemize{
#'   \item{\code{\link{fDiss}}} 
#'   \item{\code{\link{corDiss}}} 
#'   \item{\code{\link{sid}}} 
#'   \item{\code{\link{orthoDiss}}} 
#'   }
#' The functions available for evaluating similarity/dissimilarity 
#' matrices are:
#' \itemize{
#'   \item{\code{\link{simEval}}} 
#'   }
#' The functions available for removing irrelevant spectra from a 
#' reference set are:
#' \itemize{
#'   \item{\code{\link{neigCleaning}}} 
#'   }
#' The functions available for modeling spectral data using 
#' memory-based learning are:
#' \itemize{
#'   \item{\code{\link{mblControl}}}
#'   \item{\code{\link{mbl}}} 
#'   }
#' Other supplementary functions are:
#' \itemize{
#'   \item{\code{\link{print.localOrthoDiss}}}
#'   \item{\code{\link{print.mbl}}}
#'   \item{\code{\link{plot.mbl}}}
#'   \item{\code{\link{plot.orthoProjection}}}
#'   \item{\code{\link{print.orthoProjection}}}
#'   }
#' @docType package
#' @name resemble-package
#' @aliases resemble-package resemble
#' @title Overview of the functions in the resemble package
#' @import Rcpp RcppArmadillo foreach iterators
#' @useDynLib resemble
#' @author Leonardo Ramirez-Lopez \email{ramirez.lopez.leo@@gmail.com} & Antoine Stevens
######################################################################
# resemble
# Copyrigth (C) 2014 Leonardo Ramirez-Lopez and Antoine Stevens
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
######################################################################


## History:
## 09.03.2014 Leo     History comments were added to the function files

NULL
