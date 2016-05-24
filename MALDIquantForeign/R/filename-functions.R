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

.cleanFilename <- function(x) {
  gsub(pattern="([[:punct:]]|[[:space:]])+", replacement="_", x=x)
}

#' Determine file extension
#'
#' @param x \code{character}, filename.
#'
#' @return \code{character}, file extension.
#'
#' @seealso \code{\link[MALDIquant]{MassSpectrum-class}}
#' @keywords internal
#' @noRd
#' @examples
#' library("MALDIquantForeign")
#' files <- c("/home/foo/bar.txt", "foobar.pdf")
#' MALDIquantForeign:::.fileExtension(files)
#'
.fileExtension <- function(x) {
  pos <- regexpr(pattern="(\\.tar)?\\.([[:alnum:]]+)$|(/|\\\\)+[^.\\\\/]+$",
                 text=x)
  ifelse(pos > -1L, substring(x, pos+1L), x)
}

.withoutFileExtension <- function(x) {
  sub(pattern="\\.[[:alnum:]]+?$|(/|\\\\)+[^.\\\\/]+$",
      replacement="", x=.withoutCompressionExtension(x))
}

.withoutCompressionExtension <- function(x) {
  sub(pattern="\\.(zip|gz|bz2|bzip2|xz|lzma)+$", replacement="", x=x)
}

.changeFileExtension <- function(x, newExtension) {
  paste(.withoutFileExtension(x), newExtension, sep=".")
}

.cutFilenames <- function(x) {
  l <- strsplit(x, split=.Platform$file.sep, fixed=TRUE)

  nCol <- unlist(lapply(l, length))
  mCol <- max(nCol)

  m <- matrix(NA, nrow=length(x), ncol=mCol)

  for (i in seq(along=l)) {
    cols <- 1:nCol[i]
    m[i, cols] <- l[[i]]
  }

  isIdentical <- apply(m, 2, function(co)all(co[1] == co))
  isIdentical[is.na(isIdentical)] <- FALSE

  m <- as.matrix(m[, !isIdentical])

  if (length(m)) {
    filenames <- apply(m, 1, function(r) {
      do.call(file.path, as.list(na.omit(r)))
    })
  } else {
    filenames <- basename(x)
  }

  filenames
}

.uniqueBaseFilenames <- function(x, fileExtension="csv",
                                 sep="_") {
  filenames <- .cutFilenames(.withoutFileExtension(x))
  filenames <- .cleanFilename(filenames)

  empty <- nchar(filenames) <= 0
  filenames[empty] <- seq_along(empty)

  filenames <- .make.unique(filenames, sep=sep)
  paste(filenames, fileExtension, sep=".")
}

## let make unique start by 1
.make.unique <- function(x, sep="_") {
  tmp <- lapply(split(x, x), function(y) {
    n <- length(y)
    if (n > 1) {
      fmt <- paste0("%s%s%0", floor(log10(n))+1, "d")
      y <- sprintf(fmt=fmt, y, sep, 1:n)
    }
    return(y)
  })
  unsplit(tmp, x)
}
