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

#' @author Pietro Franceschi \email{pietro.franceschi@@fmach.it}, Sebastian Gibb
#' \email{mail@@sebastiangibb.de}
#' @keywords internal
#' @noRd

## original code written by Pietro Franceschi
## modified for MALDIquantForeign by Sebastian Gibb
.importCdf <- function(file, centroided=FALSE, massRange=c(0, Inf),
                          minIntensity=0, verbose=FALSE) {
  if (!requireNamespace("RNetCDF", quietly=TRUE)) {
        stop("For netCDF support the ", sQuote("RNetCDF"),
             " package is needed. \n",
             "Please install ", sQuote("RNetCDF"),
             ": install.packages(\"RNetCDF}\")")
  }

  nc <- RNetCDF::open.nc(file)
  on.exit(RNetCDF::close.nc(nc))

  scanIndex <- as.integer(RNetCDF::var.get.nc(nc, variable="scan_index"))
  scanLength <- as.integer(RNetCDF::var.get.nc(nc, variable="point_count"))
  retentionTime <- as.double(RNetCDF::var.get.nc(nc,
                               variable="scan_acquisition_time"))

  ## read and process the spectra/peaks
  l <- vector(mode="list", length=length(scanIndex))

  for (i in seq(along=l)) {
    mass <- as.double(RNetCDF::var.get.nc(nc, variable="mass_values",
                                          start=scanIndex[i]+1L,
                                          count=scanLength[i]))
    intensity <- as.double(RNetCDF::var.get.nc(nc, variable="intensity_values",
                                               start=scanIndex[i]+1L,
                                               count=scanLength[i]))
    l[[i]] <- .createMassObject(mass=mass, intensity=intensity,
                                metaData=list(file=file,
                                              number=i,
                                              retentionTime=retentionTime[i],
                                              scanIndex=scanIndex[i]),
                                centroided=centroided, massRange=massRange,
                                minIntensity=minIntensity, verbose=verbose)
  }
  l
}
