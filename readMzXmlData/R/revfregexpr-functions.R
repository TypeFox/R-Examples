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

#' Pattern matching.
#'
#' This function looks for the last match to argument \code{pattern} in a file
#' \code{file}.
#'
#' @param pattern \code{character}, string containing a regular expression.
#' @param filename \code{character}, name of file.
#'
#' @return \code{double}, position of match.
#'
#' @rdname revfregexpr 
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @keywords internal
#'
.revfregexpr <- function(pattern, file) {
  bufferSize <- 1024L
  size <- file.info(file)$size
  pos <- double()
  f <- file(file, "rb")

  if (size < bufferSize) {
    bufferSize <- size
  }

  readPos <- c(seq(from=size-bufferSize, to=0, by=-bufferSize), 0)
  n <- length(readPos)
  readBufferSize <- c(rep(bufferSize, n-1), readPos[n-1])

  for (i in seq(along=readPos)) {
    seek(f, where=readPos[i])
    
    p <- gregexpr(pattern=pattern,
                  text=readChar(f, nchars=readBufferSize[i]))[[1]]

    if (p[1] > 0) {
      pos <- readPos[i]+tail(p, 1)
      break 
    }
  }
  close(f)

  return(pos)
}

