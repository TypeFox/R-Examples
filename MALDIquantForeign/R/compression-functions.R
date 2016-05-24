## Copyright 2013 Sebastian Gibb
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

.isCompressed <- function(x) {
  z <- c("bz2", "bzip2", "gz", "lzma", "xz")
  .fileExtension(x) %in% c("zip", z, paste("tar", z, sep="."))
}

.isTar <- function(x)grepl(pattern="^tar", x=.fileExtension(x))

.isZip <- function(x)grepl(pattern="^zip$", x=.fileExtension(x))

.isPackedOrCompressed <- function(x) {
  .isCompressed(x) | .isTar(x)
}

# unpack/uncompress files and return temporary filenames
.uncompress <- function(x, verbose=TRUE) {
  f <- lapply(x, function(path) {
    if (!.isPackedOrCompressed(path)) {
      return(x)
    } else {
      if (.isTar(path)) {
        return(.unpacking(x, fun=untar, verbose=verbose))
      } else if (.isZip(path)) {
        return(.unpacking(x, fun=unzip, verbose=verbose))
      } else {
        return(.gunzip(x, verbose=verbose))
      }
    }
  })
  unlist(f)
}

# unpack and return tmp filename
.unpacking <- function(filename, destdir, fun, verbose=FALSE, ...) {
  if (missing(destdir)) {
    pattern <- paste0(.withoutFileExtension(basename(filename)), "_")
    destdir <- file.path(tempdir(), "MALDIquantForeign_uncompress",
                         tempfile(pattern=pattern, tmpdir=""))
  }

  funName <- deparse(substitute(fun))
  fun <- as.function(fun)

  .msg(verbose, funName, " ", filename, " to ", destdir, ".")

  unpacked <- fun(filename, exdir=destdir, ...)
  if (!length(unpacked)) {
    stop(funName, " failed!")
  }
  destdir
}

# gunzip and return tmp filename
.gunzip <- function(filename, destfile, verbose=FALSE) {
  if (!file.exists(filename)) {
    stop(sQuote(filename), " doesn't exist!")
  }

  if (missing(destfile)) {
    tmpdir <- file.path(tempdir(), "MALDIquantForeign_uncompress")
    if (!file.exists(tmpdir)) {
      dir.create(tmpdir, showWarnings=FALSE, recursive=TRUE)
    }

    pattern <- paste0(.withoutFileExtension(basename(filename)), "_")
    fileext <- paste0(".",
                      .fileExtension(.withoutCompressionExtension(filename)))
    destfile <- tempfile(pattern=pattern, tmpdir=tmpdir, fileext=fileext)
  }

  .msg(verbose, "gunzip ", filename, " to ", destfile, ".")

  fi <- gzfile(filename, open="rb")
  on.exit(close(fi))

  fo <- file(destfile, open="wb")
  on.exit(close(fo), add=TRUE)

  repeat {
    b <- readBin(fi, what=raw(), n=1e6) ## n==1e6 => nearly 50Mb

    if (length(b)) {
      writeBin(b, con=fo)
    } else {
      break
    }
  }

  destfile
}

.cleanupUncompressedTmpFiles <- function() {
  unlink(file.path(tempdir(), "MALDIquantForeign_uncompress"), recursive=TRUE)
}
