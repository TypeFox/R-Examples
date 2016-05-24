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

#' Parse msd files.
#'
#' This function parses msd files.
#'
#' @param file \code{character}, path to msd file.
#' @param verbose \code{logical}, verbose output?
#'
#' @return Returns a list with metadata and spectra.
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @keywords internal
#' @noRd
.parseMsd <- function(file, verbose=FALSE, ...) {
  XML::xmlEventParse(file=file,
                     handlers=.msdHandlers(fileName=file,
                                           verbose=verbose),
                     addContext=FALSE, useTagName=TRUE, useDotNames=TRUE,
                     ...)$getData()
}

#' Parse msd files.
#'
#' This function is defines handlers for XML SAX parser. Internal use only.
#'
#' @param fileName \code{character}, path to msd file
#' @param verbose \code{logical}, verbose output?
#'
#' @return function closure
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @references mMass homepage: \url{http://mmass.org/}
#' @keywords internal
#' @noRd
.msdHandlers <- function(fileName, verbose=FALSE) {
  ## define local variables

  ## handle different mzXML versions
  msdVersion <- 0

  ## save last opened tag (needed for .text()-processing)
  text <- character()
  endian <- "little"
  precision <- 32

  ## build final list
  xml <- list()
  xml$metaData <- list()
  xml$spectrum <- list()

  ## handlers for specific tags
  ## mSD
  mSD <- function(name, attrs) {
    ## fetch version information
    msdVersion <<- readMzXmlData:::.grepDouble(attrs["version"])

    .msg(verbose, "Found mSD document (version: ", msdVersion, ").")
  }

  ## mSD/description/date
  date <- function(name, attrs) {
      xml$metaData[["acquisitionDate"]] <<-
        readMzXmlData:::.attributeToString(attrs, "value")
  }

  operator <- function(name, attrs) {
      xml$metaData[["owner"]] <<-
        readMzXmlData:::.attributeToString(attrs, "value")
  }

  institution <- function(name, attrs) {
      xml$metaData[[name]] <<-
        readMzXmlData:::.attributeToString(attrs, "value")
  }

  ## mSD/spectrum
  spectrum <- function(name, attrs) {
    attributeNames <- c("scanNumber", "msLevel", "retentionTime",
                        "precursorMZ", "precursorCharge", "polarity")
    attributeNames <- intersect(attributeNames, names(attrs))
    for (i in attributeNames) {
      xml$metaData[[i]] <<- readMzXmlData:::.attributeToDouble(attrs, i)
    }
  }

  ## mSD/spectrum/mzArray
  mzArray <- function(name, attrs) {
    precision <<- readMzXmlData:::.attributeToDouble(attrs, "precision")
    endian <<- readMzXmlData:::.attributeToString(attrs, "endian")
  }

  ## default functions to catch tags without a handler
  .endElement <- function(name, attrs) {
    if (name == "title") {
      xml$metaData[c("name", "fullName", "sampleName")] <<- text
    } else if (name == "notes") {
      xml$metaData[["comment"]] <<- text
    } else if (name == "mzArray") {
      xml$spectrum[["mass"]] <<- .decodeArray()
    } else if (name == "intArray") {
      xml$spectrum[["intensity"]] <<- .decodeArray()
    }

    text <<- character()
  }

  .text <- function(x) {
    text <<- paste0(text, x)
  }

  .decodeArray <- function() {
    readMzXmlData:::.base64decode(x=text, endian=endian,
                                  size=round(precision/8L),
                                  compressionType="gzip")
  }

  ## return statement (please call getData())
  return(list(getData=function() {return(xml)},
              mSD=mSD,
              date=date,
              institution=institution,
              instrument=institution,
              spectrum=spectrum,
              mzArray=mzArray,
              intArray=mzArray,
              .endElement=.endElement,
              .text=.text))
}
