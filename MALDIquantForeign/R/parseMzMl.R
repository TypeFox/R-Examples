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

#' Parse mzML files.
#'
#' This function parses mzML files.
#'
#' @param file \code{character}, path to mzML file
#' @param verbose \code{logical}, verbose output?
#'
#' @return Returns a list with metadata and spectra.
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @keywords internal
#' @noRd
.parseMzMl <- function(file, verbose=FALSE, ...) {
  XML::xmlEventParse(file=file,
                     handlers=.mzMlHandlers(fileName=file,
                                            verbose=verbose),
                     addContext=FALSE, useTagName=TRUE, useDotNames=TRUE,
                     ...)$getData()
}

#' Parse mzML files.
#'
#' This function is defines handlers for XML SAX parser. Internal use only.
#'
#' TODO: fetch obo map file and use it
#'
#' @param fileName \code{character}, path to mzML file
#' @param verbose \code{logical}, verbose output?
#'
#' @return function closure
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @references
#' Definition of \code{mzML} format:
#' \url{http://www.psidev.info/mzml_1_0_0}
#' @keywords internal
#' @noRd
.mzMlHandlers <- function(fileName, verbose=FALSE) {
  ## define local variables

  ## handle different mzXML versions
  mzMlVersion <- 0

  ## save last opened tag (needed for .text()-processing)
  openTag <- ""

  ## store current scan values
  nSpectra <- 0L
  curSpecIdx <- 0L
  curRefId <- character()
  currentArrayContent <- character()
  fileCheckSum <- character()

  ## supported?
  supported <- TRUE

  ## build final list
  xml <- list()
  xml$metaData <- list()
  xml$spectra <- list()
  ## imzML information
  xml$ims <- list()

  ## reference values
  references <- new.env(parent=emptyenv())

  ## handlers for specific tags
  ## mzML
  mzML <- function(name, attrs) {
    ## fetch version information
    mzMlVersion <<- readMzXmlData:::.grepDouble(attrs["xsi:schemaLocation"])

    .msg(verbose, "Found mzML document (version: ", mzMlVersion, ").")
  }

  ## mzML/referenceableParamGroupList/referenceableParamGroup
  referenceableParamGroup <- function(name, attrs) {
    openTag <<- name
    curRefId <<- readMzXmlData:::.attributeToString(attrs, "id", required=TRUE)
  }

  referenceableParamGroupRef <- function(name, attrs) {
    if (!is.null(references[[attrs["ref"]]])) {
      xml$spectra[[curSpecIdx]]$metaData <<-
        modifyList(xml$spectra[[curSpecIdx]]$metaData,
                   references[[attrs["ref"]]])
    }
  }

  ## mzML/run/spectrumList
  spectrumList <- function(name, attrs) {
    ## fetch number of spectra
    nSpectra <<- readMzXmlData:::.attributeToDouble(attrs, "count",
                                                    required=TRUE)

    .msg(verbose, "Found ", nSpectra, " ",
         ifelse(nSpectra == 1, "spectrum", "spectra"), ".")
  }

  ## mzML/run/spectrumList/spectrum
  spectrum <- function(name, attrs) {
    openTag <<- name
    curSpecIdx <<- curSpecIdx + 1

    supported <<- TRUE

    xml$spectra[[curSpecIdx]] <<- list()
    xml$spectra[[curSpecIdx]]$metaData <<- list()

    xml$spectra[[curSpecIdx]]$metaData[["id"]] <<-
      readMzXmlData:::.attributeToString(attrs, "id", required=TRUE)

    xml$spectra[[curSpecIdx]]$metaData[["numberInFile"]] <<-
      curSpecIdx

    ## IMS extension
    if (length(xml$ims)) {
      xml$ims$ibd[[curSpecIdx]] <<- matrix(NA, nrow=2, ncol=3,
                                           dimnames=list(c("mass", "intensity"),
                                                         c("offset", "length",
                                                           "encodedLength")))
      xml$spectra[[curSpecIdx]]$metaData$imaging <<- list()
    }
    .msg(verbose, "Processing spectrum ", curSpecIdx, "/", nSpectra,
        " (id: ", attrs["id"], ") ...")
  }

  chromatogram <- function(name, attrs) {
    warning("<chromatogram> tag is not supported!")
    supported <<- FALSE
  }

  ## *cvParam
  cvParam <- function(name, attrs) {
    ## polarity
    if (.isAttrSet(attrs, "MS:1000129", "negative scan")) {
      .setCvValue("negative", "polarity")
      return()
    } else if (.isAttrSet(attrs, "MS:1000130", "positive scan")) {
      .setCvValue("positive", "polarity")
      return()
    }

    ## ms level
    if (.isAttrSet(attrs, "MS:1000511",  "ms level")) {
      .setCvValue(
        readMzXmlData:::.attributeToDouble(attrs, "value", required=TRUE),
        "msLevel")
      return()
    }

    ## centroid data
    if (.isAttrSet(attrs, "MS:1000127", "centroid spectrum")) {
      .setCvValue(1, "centroided")
      return()
    }

    ## precision
    if (.isAttrSet(attrs, "MS:1000521", "32-bit float")) {
      .setCvValue(32, "precision")
      return()
    } else if (.isAttrSet(attrs, "MS:1000523", "64-bit float")) {
      .setCvValue(64, "precision")
      return()
    }

    ## compression
    if (.isAttrSet(attrs, "MS:1000576", "no compression")) {
      .setCvValue("none", "compressionType")
      return()
    } else if (.isAttrSet(attrs, "MS:1000574", "zlib compression")) {
      .setCvValue("gzip", "compressionType")
      return()
    }

    ## data arrays
    if (.isAttrSet(attrs, "MS:1000514", "m/z array")) {
      .setCvValue("mass", "currentArray")
      return()
    } else if (.isAttrSet(attrs, "MS:1000515", "intensity array")) {
      .setCvValue("intensity", "currentArray")
      return()
    } else if (.isAttrSet(attrs, "MS:1000517", "signal to noise array")) {
      .setCvValue("snr", "currentArray")
      return()
    }

    ## IMS extensions
    if (.isAttrSet(attrs, "IMS:1000080", "universally unique identifier")) {
      xml$ims$uuid <<-
        gsub(pattern="[[:punct:]]", replacement="",
             x=readMzXmlData:::.attributeToString(attrs, "value",
                                                  required=TRUE))
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000090", "ibd MD5")) {
      xml$ims$md5 <<-
        readMzXmlData:::.attributeToString(attrs, "value")
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000091", "ibd SHA-1")) {
      xml$ims$sha1 <<-
        readMzXmlData:::.attributeToString(attrs, "value")
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000030", "continuous")) {
      xml$ims$type <<- "continuous"
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000031", "processed")) {
      xml$ims$type <<- "processed"
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000042", "max count of pixel x")) {
      xml$metaData$imaging$size["x"] <<-
        readMzXmlData:::.attributeToDouble(attrs, "value", required=TRUE)
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000043", "max count of pixel y")) {
      xml$metaData$imaging$size["y"] <<-
        readMzXmlData:::.attributeToDouble(attrs, "value", required=TRUE)
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000044", "max dimension x")) {
      xml$metaData$imaging$dim["x"] <<-
        readMzXmlData:::.attributeToDouble(attrs, "value", required=TRUE)
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000045", "max dimension y")) {
      xml$metaData$imaging$dim["y"]<<-
        readMzXmlData:::.attributeToDouble(attrs, "value", required=TRUE)
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000046", "pixel size x")) {
      xml$metaData$imaging$pixelSize["x"] <<-
        readMzXmlData:::.attributeToDouble(attrs, "value", required=TRUE)
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000047", "pixel size y")) {
      xml$metaData$imaging$pixelSize["y"] <<-
        readMzXmlData:::.attributeToDouble(attrs, "value", required=TRUE)
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000050", "position x")) {
      xml$spectra[[curSpecIdx]]$metaData$imaging$pos["x"] <<-
        readMzXmlData:::.attributeToDouble(attrs, "value", required=TRUE)
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000051", "position y")) {
      xml$spectra[[curSpecIdx]]$metaData$imaging$pos["y"] <<-
        readMzXmlData:::.attributeToDouble(attrs, "value", required=TRUE)
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000102", "external offset")) {
      xml$ims$ibd[[curSpecIdx]][xml$spectra[[curSpecIdx]]$metaData$currentArray,
                                "offset"] <<-
        readMzXmlData:::.attributeToDouble(attrs, "value", required=TRUE)
      return()
    }

    if (.isAttrSet(attrs, "IMS:1000103", "external array length")) {
      xml$ims$ibd[[curSpecIdx]][xml$spectra[[curSpecIdx]]$metaData$currentArray,
                                "length"] <<-
        readMzXmlData:::.attributeToDouble(attrs, "value", required=TRUE)
      return()
    }
    if (.isAttrSet(attrs, "IMS:1000104", "external encoded length")) {
      xml$ims$ibd[[curSpecIdx]][xml$spectra[[curSpecIdx]]$metaData$currentArray,
                                "encodedLength"] <<-
        readMzXmlData:::.attributeToDouble(attrs, "value", required=TRUE)
      return()
    }
  }

  ## default functions to catch tags without a handler
  .startElement <- function(name, attrs) {
    openTag <<- name
  }

  .endElement <- function(name, attrs) {
    if (name == "binary" && supported) {
      .decodeArray()
    } else if (name == "fileChecksum") {
      .calculateFileChecksum()
    }

    if (openTag == name) {
      openTag <<- ""
    }
  }

  .text <- function(x) {
    if (openTag == "binary" && supported) {
      currentArrayContent <<- paste0(currentArrayContent, x)
    } else if (openTag == "fileChecksum") {
      fileCheckSum <<- paste0(fileCheckSum, x)
    }
  }

  ## helper functions
  .isAttrSet <- function(attrs, id, name) {
    return(attrs["accession"] == id || attrs["name"] == name)
  }

  .setCvValue <- function(x, name) {
    if (openTag == "referenceableParamGroup") {
      if (is.null(references[[curRefId]])) {
        references[[curRefId]] <<- list()
      }

      references[[curRefId]][[name]] <<- x
    } else if (openTag == "spectrum" || openTag == "binaryDataArray") {
      xml$spectra[[curSpecIdx]]$metaData[[name]] <<- x
    } else {
      xml$metaData[[name]] <<- x
    }
  }

  .decodeArray <- function() {
    if (length(currentArrayContent)) {
      ## read base64 encoded array content (endian must be "little")
      content <- readMzXmlData:::.base64decode(x=currentArrayContent,
        endian="little", size=round(xml$spectra[[curSpecIdx]]$metaData$precision/8),
        compressionType=xml$spectra[[curSpecIdx]]$metaData$compressionType)

      xml$spectra[[curSpecIdx]][[xml$spectra[[curSpecIdx]]$metaData$currentArray]] <<-
        content

      ## clear array content
      currentArrayContent <<- character()
    }

    ## clear metaData
    xml$spectra[[curSpecIdx]]$metaData[c("precision", "compressionType",
                                         "currentArray")] <<- NULL
  }

  .calculateFileChecksum <- function() {
    n <- nchar(fileCheckSum)

    if (n <= 0) {
      return()
    }

    ## sha1 sum for this file (from the beginning of the file up to (and
    ## including) the opening tag <fileChecksum>
    .msg(verbose, "Look for '<fileChecksum>' position ...")

    ## 14 == nchar("<fileChecksum>")
    checkSumPos <- readMzXmlData:::.revfregexpr("<fileChecksum>", fileName) + 14

    .testChecksum(fileName, fileCheckSum, algo="sha1", length=checkSumPos-1,
                  verbose=verbose)
  }

  ## return statement (please call getData())
  return(list(getData=function() {return(xml)},
              mzML=mzML,
              referenceableParamGroup=referenceableParamGroup,
              referenceableParamGroupRef=referenceableParamGroupRef,
              spectrumList=spectrumList,
              spectrum=spectrum,
              chromatogram=chromatogram,
              cvParam=cvParam,
              .startElement=.startElement,
              .endElement=.endElement,
              .text=.text))
}
