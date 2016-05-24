## Copyright 2015 Sebastian Gibb
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

.writeImzMlDocument <- function(x, file,
                                id=.withoutFileExtension(basename(file)),
                                processed=TRUE, uuid=.uuid(),
                                coordinates=NULL, pixelSize=c(100, 100), ...) {
  if(isMassSpectrum(x)) {
    x <- list(x)
  }

  if (is.null(metaData(x[[1L]])$imaging) && is.null(coordinates)) {
    stop("The spectra contain no imaging information.")
  }

  isCoordinatesMatrix <- is.matrix(coordinates) &&
                         ncol(coordinates) == 2L &&
                         nrow(coordinates) == length(x)

  if (!is.null(coordinates) && !isCoordinatesMatrix) {
    stop("The ", sQuote("coordinates"),
         " argument has to be a matrix with two columns (x and y position)!")
  }

  if (isCoordinatesMatrix) {
    size <- apply(coordinates, 2, max)
    dimension <- size * pixelSize

    x[[1L]]@metaData$imaging <- list(size=size,
                                     dim=dimension,
                                     pixelSize=pixelSize)
    coordinates(x) <- coordinates
  }

  ibdFile <- .changeFileExtension(file, "ibd")
  .writeIbd(x, ibdFile, uuid=uuid, processed=processed)
  sha1 <- digest::digest(ibdFile, algo="sha1", file=TRUE)

  .writeMzMlDocument(x=.addIbdOffsets(x), file=file, id=id,
                     imsArgs=list(uuid=uuid, sha1=sha1, processed=processed),
                     ...)
}
