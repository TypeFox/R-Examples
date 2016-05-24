## Copyright 2013-2015 Sebastian Gibb
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

.importImzMl <- function(file, centroided=FALSE, massRange=c(0, Inf),
                         minIntensity=0, coordinates=NULL,
                         verbose=FALSE) {

  .msg(verbose, "Reading spectrum from ", sQuote(file), " ...")

  if (!file.exists(file)) {
    stop("File ", sQuote(file), " doesn't exists!")
  }

  ibdFilename <- .changeFileExtension(file, "ibd")

  if (!file.exists(ibdFilename)) {
    stop("File ", sQuote(ibdFilename), " doesn't exists!")
  }

  s <- .parseMzMl(file=file, verbose=verbose)

  ibd <- file(ibdFilename, open="rb")
  on.exit(close(ibd))

  ## test UUID
  uuid <- paste0(readBin(ibd, raw(), n=16, size=1, signed=TRUE, endian=
                         "little"), collapse="")

  if (is.null(s$ims$uuid)) {
    warning("There is not any UUID in ", sQuote(file), "!")
  } else if (tolower(uuid) != tolower(s$ims$uuid)) {
    warning("The UUID in ", sQuote(file) , " and ", sQuote(ibdFilename),
            "do not match!")
  } else if (!.isUuidV4(uuid)) {
    warning("The UUID: ", uuid, " is not valid!")
  }

  ## test checksums
  if (!is.null(s$ims$md5)) {
    .testChecksum(ibdFilename, s$ims$md5, algo="md5", verbose=verbose)
  } else if (!is.null(s$ims$sha1)) {
    .testChecksum(ibdFilename, s$ims$sha1, algo="sha1", verbose=verbose)
  } else {
    warning("At least one checksum (SHA-1 or MD5) for the idb file ",
            "should be provided in the imzML file.")
  }

  sel <- seq_along(s$ims$ibd)

  isCoordinatesMatrix <- is.matrix(coordinates) && ncol(coordinates) == 2L

  if (!is.null(coordinates) && !isCoordinatesMatrix) {
    stop("The ", sQuote("coordinates"),
         " argument has to be a matrix with two columns (x and y position)!")
  }

  if (isCoordinatesMatrix) {
    pos <- do.call(rbind, lapply(s$spectra, function(x)x$metaData$imaging$pos))
    sel <- match(paste(coordinates[, 1L], coordinates[, 2L], sep=":"),
                 paste(pos[, 1L], pos[, 2L], sep=":"))
    if (anyNA(sel)) {
      warning("The following rows contain invalid coordinates: ",
              paste(which(is.na(sel)), collapse=", "))

      sel <- sort(sel[!is.na(sel)])
    }
  }

  .readValues <- function(file, x, column, isSeekNeeded) {
    if (isSeekNeeded) {
      ## WARNING: we know that `seek` is discouraged on some platforms,
      ## namely Windows. See `?seek` for details.
      seek(file, where=x[column, "offset"])
    }
    n <- x[column, "length"]
    e <- x[column, "encodedLength"]
    return(readBin(file, double(), n=n, size=e/n, signed=TRUE, endian="little"))
  }

  n <- length(sel)
  spectra <- vector(mode="list", length=n)

  isProcessed <- s$ims$type == "processed"
  isSeekNeeded <- length(s$ims$ibd) > length(sel)

  if (!isProcessed) {
    mass <- .readValues(ibd, s$ims$ibd[[sel[1L]]], "mass", isSeekNeeded)
  }

  ## read mass and intensity values
  for (i in seq(along=sel)) {
    .msg(verbose, "Reading binary data for spectrum ", i, "/", n, " ...")

    m <- modifyList(s$metaData, s$spectra[[sel[i]]]$metaData)
    m$file <- file

    if (isProcessed) {
      mass <- .readValues(ibd, s$ims$ibd[[sel[i]]], "mass", isSeekNeeded)
    }
    intensity <- .readValues(ibd, s$ims$ibd[[sel[i]]], "intensity", isSeekNeeded)
    spectra[[i]] <- .createMassObject(mass=mass, intensity=intensity,
                                      metaData=m, centroided=centroided,
                                      massRange=massRange,
                                      minIntensity=minIntensity,
                                      verbose=verbose)
  }
  spectra
}
