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

#' Converts double to single precision.
#'
#' This function simulates the conversion of floating point numbers from double
#' precision (64bit, R: \code{double()}, C: \code{double}) to single precision
#' (32bit, R: \code{none}, C: \code{float}). It follows IEEE 754 standard.
#'
#' @param x \code{double}, numeric value which should reduce from double
#'  precision to single precision.
#' @return \code{double} (in single precision).
#'
#' @details
#' The same could be done in C by using casts:
#' \preformatted{
#' double precision32(double value) {
#'   float x=value;
#'   return (double)x;
#' }}
#' @seealso \code{\link[readBrukerFlexData]{.changePrecision}},
#'  \code{\link[readBrukerFlexData]{.hpc}}
#' @rdname double2singlePrecision
#' @keywords internal
#' @references
#' IEEE 754 standard: \url{http://754r.ucbtest.org/standards/754.pdf}
#' @examples
#' ## load library
#' library("readBrukerFlexData")
#'
#' ## show more details
#' oldDigits <- options()$digits
#' options(digits=22)
#'
#' ## a test number
#' num <- 1/3
#'
#' num
#' readBrukerFlexData:::.double2singlePrecision(num)
#'
#' ## reset digits option
#' options(digits=oldDigits)
#'
.double2singlePrecision <- function(x) {
  stopifnot(is.double(x))
  return(.changePrecision(x, size=4))
}

#' Change precision.
#'
#' This function converts double values to double values in a given
#' precision (only correctly working for cut a higher precision to a lower one;
#' e.g.  IEEE 754 double precision to IEEE 754 single precision)
#'
#' @param x \code{double}, a vector of double values which should converted.
#' @param size \code{integer}, how many bytes using for target precision
#'  e.g. IEEE 754 single precision: \code{size=4},
#'  IEEE 754 double precision: \code{size=8}
#' @return \code{double}.
#' @seealso \code{\link[readBrukerFlexData]{.double2singlePrecision}},
#' @rdname changePrecision
#' @keywords internal
#'
.changePrecision <- function(x, size) {
  ## create a raw object to avoid direct file access
  virtualCon <- raw()
  ## write binary data to raw object and change (mostly cut) precision to size
  ## size==4 # 32bit, single precision
  ## size==8 # 64bit, double precision
  virtualCon <- writeBin(object=x, con=virtualCon, size=size)
  ## re-read data
  x <- readBin(con=virtualCon, what=double(), size=size, n=length(x))
  return(x)
}

