## Copyright 2013-2014 Sebastian Gibb
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

.importAnalyze <- function(file, centroided=FALSE, massRange=c(0, Inf),
                           minIntensity=0, verbose=FALSE) {
  header <- .readAnalyzeHdr(file, verbose=verbose)
  mass <- .readAnalyzeMass(.changeFileExtension(file, "t2m"),
                           header=header,
                           verbose=verbose)

  massRange <- MALDIquant:::.reorderRange(massRange)
  massIdx <- which(massRange[1] <= mass & mass <= massRange[2])
  skip <- c(massIdx[1]-1, length(mass)-massIdx[length(massIdx)])
  mass <- mass[massIdx]

  intensity <- .readAnalyzeIntensity(.changeFileExtension(file, "img"),
                                     header=header,
                                     ni=length(mass),
                                     skip=skip,
                                     verbose=verbose)

  l <- vector(mode="list", length=header$nx*header$ny)

  for (x in 1:header$nx) {
    for (y in 1:header$ny) {
      l[[(x-1)*header$ny+y]] <-
        .createMassObject(mass=mass, intensity=intensity[, x, y],
                          metaData=list(file=file,
                                        imaging=list(pos=c(x=x, y=y),
                                                     pixelSize=c(x=header$xd,
                                                                  y=header$yd))),
                          centroided=centroided, massRange=massRange,
                          minIntensity=minIntensity, verbose=verbose)
    }
  }
  l
}
