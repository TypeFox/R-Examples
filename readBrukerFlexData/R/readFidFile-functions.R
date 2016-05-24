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

#' Reads binary fid file.
#'
#' This function reads a binary fid file. A fid file contains intensities for
#' all measured time points.
#'
#'
#' @param fidFile \code{character} path to fid file e.g.
#'  Pankreas_HB_L_061019_A10/0_a19/1/1SLin/fid
#' @param nIntensities number of data entries
#'  (total count; get from acqu file)
#' @param endian endianness of the fid file
#'  (\sQuote{little} or \sQuote{big}; default: \sQuote{little})
#'
#' @return
#'  A vector of intensity values.
#' @seealso
#'  \code{\link[readBrukerFlexData]{.readAcquFile}},
#'  \code{\link[readBrukerFlexData]{readBrukerFlexFile}}
#' @keywords internal IO
#' @rdname readFidFile
#'
.readFidFile <- function(fidFile, nIntensities, endian="little") {
  if (!file.exists(fidFile)) {
    stop("File ", sQuote(fidFile), " doesn't exists!")
  }

  con <- file(fidFile, "rb")
  intensity <- readBin(con, integer(), n=nIntensities, size=4, endian=endian)
  close(con)

  return(as.double(intensity))
}

