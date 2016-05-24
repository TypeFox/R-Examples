## Copyright 2014-2015 Sebastian Gibb
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

#' Create a MassSpectrum or a MassPeaks object.
#'
#' This function creates a MassSpectrum or MassPeaks object (depending on the
#' centroided argument).
#'
#' @param mass mass data
#' @param intensity intensity data
#' @param snr snr data
#' @param metaData \code{list}, metaData
#' @param centroided \code{logical}, centroided
#'  (if TRUE => MassPeaks, if FALSE => MassSpectrum)
#' @param massRange \code{double}, length == 2, trim spectrum to
#'  \code{massRange}.
#' @param minIntensity \code{double}, minimal intensity
#' @param verbose \code{logical}, verbose output?
#'
#' @return Returns a MassSpectrum or a MassPeaks object.
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @rdname createMassObject
#' @keywords internal
#' @noRd
.createMassObject <- function(mass, intensity, snr=NULL, metaData=list(),
                              centroided=FALSE,
                              massRange=c(0, Inf),
                              minIntensity=0,
                              verbose=FALSE) {

  if (!is.null(metaData$centroided) ||
      !is.null(metaData$dataProcessing$centroided)) {
    isCentroided <-
      isTRUE(as.logical(as.numeric(metaData$centroided))) |
      isTRUE(as.logical(as.numeric(metaData$dataProcessing$centroided)))

    if (isCentroided != centroided) {
      warning("According to the metadata information the imported data are ",
              ifelse(isCentroided, "", "not "), "centroided, ",
              "but they are treated as ",
              ifelse(centroided, "centroided (MassPeaks)",
                     "profile (MassSpectrum)"), " data. Maybe you want to use ",
              sQuote(paste0("centroided=", as.character(isCentroided))),
              ". See ", sQuote("?import"), " for details.", immediate.=TRUE)
    }
  }

  .msg(centroided & verbose,
       "Assume centroided data and creating a MassPeaks object.")

  ## trim AbstractMass object
  massRange <- MALDIquant:::.reorderRange(massRange)

  ## we don't use MALDIquant::trim here because we want to filter on the
  ## intensity as well
  i <- which(massRange[1L] <= mass & mass <= massRange[2L] &
             intensity >= minIntensity)

  if (is.null(snr)) {
    snr <- rep.int(NA_real_, length(i))
  }

  ## create a MassPeaks object for centroided data
  if (centroided) {
    m <- createMassPeaks(mass=mass[i], intensity=intensity[i],
                         snr=snr[i], metaData=metaData)
  } else {
    m <- createMassSpectrum(mass=mass[i], intensity=intensity[i],
                            metaData=metaData)
  }
  m
}
