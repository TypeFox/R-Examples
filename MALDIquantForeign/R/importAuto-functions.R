## Copyright 2012-2014 Sebastian Gibb
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

.importAuto <- function(path, excludePattern=NULL,
                        removeEmptySpectra=TRUE, centroided=FALSE,
                        massRange=c(0, Inf), minIntensity=0,
                        verbose=FALSE,
                        ...) {

  files <- lapply(importFormats$pattern, .files, path=path,
                  excludePattern=excludePattern)
  names(files) <- importFormats$type

  ## test xml files for ciphergen format
  files$ciphergen <- .testCiphergenXml(files$ciphergen)

  n <- vapply(files, length, integer(1))

  if (all(n)) {
    stop("Could not detect any supported file type.")
  }

  m <- which.max(n)

  .msg(verbose,
       n[m], " files of type=", sQuote(importFormats$type[m]), " found.")

  import(path=files[[m]], type=importFormats$type[m],
         pattern=importFormats$pattern[m],
         removeEmptySpectra=removeEmptySpectra, centroided=centroided,
         massRange=massRange, minIntensity=minIntensity, verbose=verbose, ...)
}

# test xml for ciphergen format
# returns files in ciphergen xml format
.testCiphergenXml <- function(files) {
  ## read first 4 lines of each file
  l <- lapply(files, readLines, n=4)
  p <- lapply(l, grepl, pattern="<spectrum>|<fileVersion>|<spectrumName>")
  s <- vapply(p, sum, integer(1))
  isCiphergen <- s >= 2
  files[isCiphergen]
}
