## Copyright 2010-2012 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This file is part of readBrukerFlexData for R and related languages.
##
## readBrukerFlexData is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## readBrukerFlexData is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with readBrukerFlexData. If not, see <http://www.gnu.org/licenses/>

#' Calculates mass from time-of-flight values.
#'
#' This function calculates mass from time of flight values
#' based on the following article:
#'  Titulaer, M. K. and Siccama, I. and Dekker, J. L. and 
#'  van Rijswijk, A. L. and Heeren R. M. and Sillevis Smitt, P. A. and
#'  Luider, T. M. (2006) \cr
#'  \dQuote{A database application for pre-processing, storage and comparison
#'  of mass spectra derived from patients and controls},
#'  \emph{BMC Bioinformatics}, \bold{7}: 403
#'  \url{http://www.ncbi.nlm.nih.gov/pubmed/16953879}
#'
#' @details Arguments are imported from metadata (acqu-file).
#' 
#' @param tof \code{double}, vector with times-of-flight 
#' @param c1 metaData$calibrationConstants[1]
#' @param c2 metaData$calibrationConstants[2]
#' @param c3 metaData$calibrationConstants[3]
#' 
#' @return \code{double}, calculated mass
#' @keywords internal
#' @rdname tof2mass
#' @references
#'  Titulaer, M. K. and Siccama, I. and Dekker, J. L. and 
#'  van Rijswijk, A. L. and Heeren R. M. and Sillevis Smitt, P. A. and
#'  Luider, T. M. (2006) \cr
#'  \dQuote{A database application for pre-processing, storage and comparison
#'  of mass spectra derived from patients and controls},
#'  \emph{BMC Bioinformatics}, \bold{7}: 403
#'  \url{http://www.ncbi.nlm.nih.gov/pubmed/16953879}
#'
.tof2mass <- function(tof, c1, c2, c3) {
  A <- c3
  B <- sqrt(1e+12/c1)
  C <- c2 - tof
  
  if (A == 0) {
    ## linear: 0 = B * sqrt(m/z) + C(times)
    return((C * C)/(B * B))
  } else {
    ## quadratic: 0 = A * (sqrt(m/z))^2 + B * sqrt(m/z) + C(times)
    return(((-B + sqrt((B * B) - (4 * A * C)))/(2 * A))^2)
  }
}
 
