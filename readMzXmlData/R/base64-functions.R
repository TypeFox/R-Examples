## Copyright 2011-2012 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This file is part of readMzXmlData for R and related languages.
##
## readMzXmlData is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## readMzXmlData is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with readMzXmlData. If not, see <http://www.gnu.org/licenses/>

#' Converts base64 character to double.
#'
#' This function converts a base64 encoded \code{character} vector to a
#' \code{double} vector.
#'
#' @param x \code{character}, base64 encoded string.
#' @param size \code{integer}, number of bytes per element in the byte stream
#'  (see \code{size} in \code{\link[base]{readBin}}).
#' @param endian \code{character}, the endian-ness
#'  (see \code{endian} in \code{\link[base]{readBin}}).
#' @param compressionType \code{character}, type of compression to use for
#'  decompression of \code{x} (see \code{type} in
#'  \code{\link[base]{memDecompress}}.
#' @return Vector of type \code{double}.
#' @rdname base64-decode
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @seealso \code{\link[base64enc]{base64decode}} from \pkg{base64enc} package
#' @keywords internal
#
.base64decode <- function(x, size, endian=.Platform$endian,
                          compressionType=c("none", "gzip")) {
  x <- base64enc::base64decode(x)
  r <- as.raw(x)

  compressionType <- match.arg(compressionType, choices=c("none", "gzip"),
                               several.ok=FALSE)
  r <- memDecompress(from=r, type=compressionType)

  return(readBin(r, what="double", n=length(r)%/%size, size=size, signed=TRUE,
                 endian=endian))
}

