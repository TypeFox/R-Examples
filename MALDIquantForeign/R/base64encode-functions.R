## Copyright 2012 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This file is part of MALDIquantForeign for R and related languages.
##
## MALDIquantForeign is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## MALDIquantForeign is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with MALDIquantForeign. If not, see <http://www.gnu.org/licenses/>

#' Converts double to base64 character.
#'
#' This function converts a \code{double} vector to a base64 encoded
#' \code{character} vector.
#'
#' @param x \code{double}, vector
#' @param size \code{integer}, number of bytes per element in the byte stream
#'  (see \code{size} in \code{\link[base]{writeBin}}).
#' @param endian \code{character}, the endian-ness
#'  (see \code{endian} in \code{\link[base]{writeBin}}).
#' @param compressionType \code{character}, type of compression to use for
#'  compression of \code{x} (see \code{type} in
#'  \code{\link[base]{memCompress}}.
#' @return Vector of type \code{character}.
#' @rdname base64-encode
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @seealso \code{\link[base64enc]{base64encode}} from \pkg{base64enc} package
#' @keywords internal
#'
.base64encode <- function(x, size, endian=.Platform$endian,
                          compressionType=c("none", "gzip")) {
  x <- writeBin(as.double(x), con=raw(), size=size, endian=endian)

  compressionType <- match.arg(compressionType, choices=c("none", "gzip"),
                               several.ok=FALSE)
  x <- memCompress(from=x, type=compressionType)

  base64enc::base64encode(x)
}
