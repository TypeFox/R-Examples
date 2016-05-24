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

#' Parse Ciphergen XML files.
#'
#' This function parses Ciphergen XML files.
#'
#' @param file \code{character}, path to Ciphergen XML file
#' @param verbose \code{logical}, verbose output?
#'
#' @return Returns a list with metadata and spectra.
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @keywords internal
#' @noRd
.parseCiphergenXml <- function(file, ...) {

  ## metaData
  listNames <- c("name",
                 "ionSourceVoltage",
                 "digitizerRate",
                 "massCalibrationA",
                 "massCalibrationB",
                 "massCalibrationT0",
                 "massCalibrationInfo",
                 "spotCorrectionFactor")
  xpath <- paste0("//spectrum/",
                  c("spectrumName",
                    paste0("acquisitionInfo/setting/", listNames[2:3]),
                    paste0("processingParameters/massCalibration/",
                           listNames[-c(1:3)])),
                 "/text()")

  doc <- xmlParse(file, ...)
  metaData <- XML::xpathApply(doc=doc, path=xpath, fun=XML::xmlValue)
  names(metaData) <- listNames
  metaData[listNames[-c(1, 7)]] <- lapply(metaData[-c(1, 7)], as.double)
  metaData$file <- normalizePath(file)

  intensity <- XML::xpathSApply(doc=doc,
                                path="//spectrum/tofData/tofDataSamples/text()",
                                fun=XML::xmlValue)
  intensity <- as.double(strsplit(intensity, "[[:space:]]+")[[1]])
  intensity <- intensity[!is.na(intensity)]

  ## tof2mass
  ## mass = U*(A*(tof-t0)^2 + B)
  ## Calibration formula was taken from:
  ## http://bioinformatics.mdanderson.org/Supplements/Datasets/KuererQC/scripts.zip
  ## file: ciphergenXMLreader.pl
  n <- 0:(length(intensity)-1)
  time <- metaData$spotCorrectionFactor*n/metaData$digitizerRate
  tof <- time-metaData$massCalibrationT0
  mass <- metaData$ionSourceVoltage *
            (sign(tof) * metaData$massCalibrationA * tof*tof +
             metaData$massCalibrationB)

  list(spectrum=list(mass=mass, intensity=intensity), metaData=metaData)
}
