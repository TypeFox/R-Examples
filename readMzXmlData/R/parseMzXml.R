## Copyright 2011-2012 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This file is part of readMzXmlData for R and related languages.
##
## readMzXmlData is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## readMzXmlData is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with readMzXmlData. If not, see <http://www.gnu.org/licenses/>

#' Parse mzXML files.
#'
#' This function parses mzXML files.
#'
#' @param file \code{character}, path to mzXML file
#' @param verbose \code{logical}, verbose output?
#'
#' @return Returns a list with metadata and spectra.
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @rdname parseMzXml
#' @keywords internal
#'
.parseMzXml <- function(file, verbose=FALSE, ...) {
  return(XML::xmlEventParse(file=file,
                            handlers=.mzXmlHandlers(fileName=file,
                                                    verbose=verbose),
                            addContext=FALSE, useTagName=TRUE, useDotNames=TRUE,
                            ...)$getData())
}

#' Parse mzXML files.
#'
#' This function is defines handlers for XML SAX parser. Internal use only.
#'
#' @param fileName \code{character}, path to mzXML file
#' @param verbose \code{logical}, verbose output?
#'
#' @return function closure
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @references
#' Definition of \code{mzXML} format:
#' \url{http://tools.proteomecenter.org/mzXMLschema.php}
#' @rdname mzXmlHandlers
#' @keywords internal
#'
.mzXmlHandlers <- function(fileName, verbose=FALSE) {
  ## define local variables

  ## handle different mzXML versions
  mzXmlVersion <- 0

  ## save last opened tag (needed for .text()-processing)
  openTag <- ""

  ## store current scan values
  nScans <- 0
  currentScanId <- 0
  currentPeaks <- ""
  precision <- 0
  byteOrder <- ""
  compressionType <- "none"

  ## sha1 tmp values
  sha1Sums <- character()
  currentSha1Id <- 0

  ## build final list
  xml <- list()
  xml$metaData <- list()
  xml$scans <- list()

  ## handlers for specific tags
  ## mzXML
  mzXML <- function(name, attrs) {
    ## fetch version information
    mzXmlVersion <<- .grepDouble(attrs["xmlns"])

    if (verbose) {
      message("Found mzXML document (version: ", mzXmlVersion, ").")
    }
  }

  ## mzXML/msRun
  msRun <- function(name, attrs) {
    ## fetch metadata for a measurement

    ## number of scans in mzXML file
    n <- .attributeToDouble(attrs, "scanCount", required=FALSE)

    nScans <<- ifelse(length(n), n, 1)

    ## optional attributes
    optAttrs <- c("startTime", "endTime")
    m <- lapply(optAttrs, .attributeTimeToDouble, attributes=attrs)
    names(m) <- optAttrs

    m <- m[!is.na(m)]

    xml$metaData <<- c(scanCount=nScans, m)
  }

  ## mzXML/msRun/parentFile
  parentFile <- function(name, attrs) {
    if (is.null(xml$metaData$parentFile)) {
      xml$metaData$parentFile <<- list()
    }

    reqAttrs <- c("fileName", "fileType", "fileSha1")
    p <- lapply(reqAttrs, .attributeToString, attributes=attrs)
    names(p) <- reqAttrs

    l <- length(xml$metaData$parentFile)+1
    xml$metaData$parentFile[[l]] <<- p
  }

  ## mzXML/msRun/msInstrument - children
  msInstrumentChild <- function(name, attrs) {
    if (is.null(xml$metaData$msInstrument)) {
      xml$metaData$msInstrument <<- list()
    }

    xml$metaData$msInstrument[[name]] <<-
      .attributeToString(attrs, "value")
  }

  ## mzXML/msRun/msInstrument/software or
  ## mzXML/msRun/dataProcessing/software
  software <- function(name, attrs) {
    reqAttrs <- c("type", "name", "version")
    optAttrs <- c("completionTime")

    if (openTag == "msInstrument") {
      xml$metaData$msInstrument$software <<- list()
      for (i in reqAttrs) {
        xml$metaData$msInstrument$software[[i]] <<-
          .attributeToString(attrs, i, required=TRUE)
      }
      for (i in optAttrs) {
        a <- .attributeToString(attrs, i)
        if (!is.na(a)) {
          xml$metaData$msInstrument$software[[i]] <<- a
        }
      }
    } else if (openTag == "dataProcessing") {
      if (is.null(xml$metaData$dataProcessing$software)) {
        xml$metaData$dataProcessing$software <<- list()
      }

      l <- length(xml$metaData$dataProcessing$software)+1
      xml$metaData$dataProcessing$software[[l]] <<- list()
      for (i in reqAttrs) {
        xml$metaData$dataProcessing$software[[l]][[i]] <<-
          .attributeToString(attrs, i, required=TRUE)
      }
      for (i in optAttrs) {
        a <- .attributeToString(attrs, i)
        if (!is.na(a)) {
          xml$metaData$dataProcessing$software[[l]][[i]] <<- a
        }
      }
    }
  }

  ## mzXML/msRun/msInstrument/operator
  operator <- function(name, attrs) {
    reqAttrs <- c("first", "last")
    optAttrs <- c("email", "phone", "URI")
    xml$metaData$msInstrument$operator <<- list()

    for (i in reqAttrs) {
      xml$metaData$msInstrument$operator[[i]] <<-
        .attributeToString(attrs, i, required=TRUE)
    }
    for (i in optAttrs) {
      a <- .attributeToString(attrs, i)
      if (!is.na(a)) {
        xml$metaData$msInstrument$operator[[i]] <<- a
      }
    }
  }

  ## mzXML/msRun/dataProcessing
  dataProcessing <- function(name, attrs) {
    if (is.null(xml$metaData$dataProcessing)) {
      xml$metaData$dataProcessing <<- list()
    }

    openTag <<- "dataProcessing"
    optAttrs <- c("intensityCutoff", "centroided", "chargeDeconvoluted",
                  "deisotoped", "spotIntegration")

    for (i in optAttrs) {
      a <- .attributeToDouble(attrs, i)
      if (length(a) && !is.na(a)) {
        xml$metaData$dataProcessing[[i]] <<- a
      }
    }
  }

  ## mzXML/msRun/dataProcessing/processingOperation
  processingOperation <- function(name, attrs) {
    optAttrs <- c("name", "value", "type")

    operations <- list()

    for (i in optAttrs) {
      a <- .attributeToString(attrs, i)
      if (!is.na(a)) {
        operations[[i]] <- a
      }
    }

    if (length(operations)) {
      if (is.null(xml$metaData$dataProcessing$operations)) {
        xml$metaData$dataProcessing$operations <<- list()
      }
      l <- length(xml$metaData$dataProcessing$operations)
      xml$metaData$dataProcessing$operations[[l+1]] <<- operations
    }
  }

  ## mzXML/msRun/scan
  scan <- function(name, attrs) {
    currentScanId <<- currentScanId + 1

    reqAttrs <- c("num", "peaksCount")
    optAttrsStr <- c("polarity", "scanType", "collisionGas", "filterLine")
    optAttrsDouble <- c("centroided", "chargeDeconvoluted", "deisotoped",
                        "ionisationEnergy", "collisionEnergy",
                        "collisionGasPressure", "cidGasPressure",
                        "startMz", "endMz", "lowMz",
                        "highMz", "basePeakMz", "basePeakIntensity",
                        "totIonCurrent", "msInstrumentID",
                        "compensationVoltage")
    optAttrsTime <- c("retentionTime")

    xml$scans[[currentScanId]] <<- list()
    xml$scans[[currentScanId]]$metaData <<- list()

    for (i in reqAttrs) {
      xml$scans[[currentScanId]]$metaData[[i]] <<-
          .attributeToDouble(attrs, i, required=TRUE)
    }

    ## msLevel is required in mzXML specification
    ## to read malformed files missing is allowed
    msLevel <- .attributeToDouble(attrs, "msLevel", required=FALSE)

    if (is.na(msLevel)) {
      msLevel <- 1
      warning("Malformed mzXML: 'msLevel' attribute is missing!\n",
              "'msLevel' == 1 is assumed.")
    }

    xml$scans[[currentScanId]]$metaData[["msLevel"]] <<- msLevel

    for (i in optAttrsStr) {
      a <- .attributeToString(attrs, i)
      if (!is.na(a)) {
        xml$scans[[currentScanId]]$metaData[[i]] <<- a
      }
    }
    for (i in optAttrsDouble) {
      a <- .attributeToDouble(attrs, i)
      if (!is.na(a)) {
        xml$scans[[currentScanId]]$metaData[[i]] <<- a
      }
    }
    for (i in optAttrsTime) {
      a <- .attributeTimeToDouble(attrs, i)
      if (!is.na(a)) {
        xml$scans[[currentScanId]]$metaData[[i]] <<- a
      }
    }

    if (verbose) {
      message("Processing scan ", currentScanId, "/", nScans,
              " (msLevel: ", attrs["msLevel"],
              ", number of peaks: ", attrs["peaksCount"], ") ...")
    }
  }

  ## mzXML/msRun/scan/scanOrigin
  scanOrigin <- function(name, attrs) {
    if (is.null(xml$scan[[currentScanId]]$metaData$scanOrigin)) {
      xml$scan[[currentScanId]]$metaData$scanOrigin <<- list()
    }

    o <- list()
    o$parentFileID <- .attributeToString(attrs, "parentFileID", required=TRUE)
    o$num <- .attributeToDouble(attrs, "num", required=TRUE)

    xml$scan[[currentScanId]]$metaData$scanOrigin[[o$num]] <<- o
  }

  ## mzXML/msRun/scan/precursorMz
  precursorMz <- function(name, attrs) {
    openTag <<- "precursorMz"

    optAttrsStr <- c("activationMethod")
    optAttrsDouble <- c("precursorCharge", "precursorIntensity",
                        "precursorScanNum", "windowWideness")
    xml$scans[[currentScanId]]$metaData$precursor <<- list()

    for (i in optAttrsStr) {
      a <- .attributeToString(attrs, i)
      if (!is.na(a)) {
        xml$scans[[currentScanId]]$metaData$precursor[[i]] <<- a
      }
    }
    for (i in optAttrsDouble) {
      a <- .attributeToDouble(attrs, i)
      if (!is.na(a)) {
        xml$scans[[currentScanId]]$metaData$precursor[[i]] <<- a
      }
    }
  }

  ## mzXML/msRun/scan/maldi
  maldi <- function(name, attrs) {
    maldi <- list()
    reqAttrs <- c("plateID", "spotID")
    optAttrs <- c("laserShootCount", "laserFrequency", "laserIntensity")

    for (i in reqAttrs) {
      maldi[[i]] <- .attributeToString(attrs, i, required=TRUE)
    }
    for (i in optAttrs) {
      a <- .attributeToString(attrs, i)
      if (!is.na(a)) {
        maldi[[i]] <- .attributeToString(attrs, i)
      }
    }

    xml$scans[[currentScanId]]$metaData$maldi <<- maldi
  }

  ## mzXML/msRun/scan/peaks
  peaks <- function(name, attrs) {
    openTag <<- "peaks"

    precision <<- .attributeToDouble(attrs, "precision", required=TRUE)
    byteOrder <<- .attributeToString(attrs, "byteOrder", required=TRUE)

    xml$scans[[currentScanId]]$metaData$precision <<- precision
    xml$scans[[currentScanId]]$metaData$byteOrder <<- byteOrder

    ## handle different versions
    if (mzXmlVersion < 3) {
      pairOrder <- .attributeToString(attrs, "pairOrder", required=TRUE)

      if (pairOrder != "m/z-int") {
        stop("Malformed mzXML: incorrect 'pairOrder' attribute ",
             "of 'peaks' field!")
      }

      xml$scans[[currentScanId]]$metaData$pairOrder <<- pairOrder

    } else {
      contentType <- .attributeToString(attrs, "contentType", required=TRUE)

      validContentTypes <- c("m/z-int", "m/z", "m/z ruler", "TOF", "intensity",
                             "S/N", "charge")

      if (!contentType %in% validContentTypes) {
        stop("Malformed mzXML: incorrect 'contentType' attribute ",
             "of 'peaks' field!")
      }

      xml$scans[[currentScanId]]$metaData$contentType <<- contentType

      ## compression
      compressionType <<- .attributeToString(attrs, "compressionType",
                                             required=TRUE)

      validCompressionTypes <- c("none", "zlib")

      if (!compressionType %in% validCompressionTypes) {
        stop("Malformed mzXML: incorrect 'compressionType' attribute ",
             "of 'peaks' field!")
      }

      compressedLen <- .attributeToDouble(attrs, "compressedLen",
                                          required=TRUE)

      xml$scans[[currentScanId]]$metaData$compressionType <<- compressionType
      xml$scans[[currentScanId]]$metaData$compressedLen <<- compressedLen

      compressionType <<- ifelse(compressionType == "zlib", "gzip", "none")
    }
  }

  ## mzXML/msRun/sha1 or
  ## mzXML/sha1
  sha1 <- function(name, attrs) {
    openTag <<- "sha1"
    currentSha1Id <<- currentSha1Id+1
    sha1Sums[currentSha1Id] <<- ""
  }

  ## default functions to catch tags without a handler
  .startElement <- function(name, attrs) {
    if (openTag == "msInstrument") {
      return(msInstrumentChild(name=name, attrs=attrs))
    } else {
      openTag <<- name
    }
  }

  .endElement <- function(name, attrs) {
    if (name == "peaks") {
      decodePeaks()
      ## clear peaks
      currentPeaks <<- ""
    }

    if (openTag == name) {
      openTag <<- ""
    }
  }

  .text <- function(x) {
    if (openTag == "precursorMz") {
      xml$scans[[currentScanId]]$metaData$precursorMz <<- as.double(x)
    } else if (openTag == "peaks") {
      currentPeaks <<- paste0(currentPeaks, x)
    } else if (openTag == "sha1") {
      sha1Sums[currentSha1Id] <<- paste0(sha1Sums[currentSha1Id], x)
    }
  }

  ## helper functions
  decodePeaks <- function() {
    ## read base64 encoded peak information
    if (byteOrder == "network") {
      endian <- "big"
    } else {
      endian <- "little"
    }

    size <- round(precision/8)

    if (size != 4 && size != 8) {
      stop("Malformed mzXML: incorrect 'precision' attribute of ",
           "'peaks' field!")
    }

    peaksCount <- xml$scans[[currentScanId]]$metaData$peaksCount

    if (peaksCount>0) {
      p <- .base64decode(x=currentPeaks, endian=endian, size=size,
                         compressionType=compressionType)
      np <- length(p) %/% 2

      if (np != peaksCount) {
        stop("Malformed mzXML: incorrect 'peakCount' attribute of ",
             "'peaks' field: expected ", peaksCount, ", found ",
             np, "  ", (3*((nchar(currentPeaks)*size)/4))/2, " (scan #",
             currentScanId, ")")
      }

      dim(p) <- c(2, np)
      mass <- p[1,]
      intensity <- p[2,]
    } else {
      mass <- double()
      intensity <- double()
    }

    xml$scans[[currentScanId]]$spectrum <<- list(mass=mass, intensity=intensity)
  }

  calculateSha1Sums <- function() {
    n <- length(sha1Sums)
    if (n <= 0) {
      return()
    }

    ## sha1 sum for this file (from the beginning of the file up to (and
    ## including) the opening tag <sha1>
    if (verbose) {
      message("Look for '<sha1>' positions ...")
    }
    sha1Pos <- .revfregexpr("<sha1>", fileName) + 6 # 6 == length("<sha1>")
    ## multiple sha1 sections are possible
    for (i in n) {
      if (verbose) {
        message("Calculate sha1-sum (", i, "/", n, ") for ", sQuote(fileName),
                ": ", appendLF=FALSE)
      }

      sha1Calc <- digest::digest(fileName, algo="sha1", file=TRUE,
                                 length=sha1Pos[i]-1)

      if (verbose) {
        message(sha1Calc)
      }

      if (sha1Sums[i] != sha1Calc) {
        warning("Stored and calculated Sha-1 sums do not match ",
                "(stored ", sQuote(sha1Sums[i]), " calculated ",
                sQuote(sha1Calc), ")!")
      }
    }
  }

  ## return statement (please call getData())
  return(list(getData=function() {return(xml)},
              mzXML=mzXML,
              msRun=msRun,
              parentFile=parentFile,
              software=software,
              operator=operator,
              dataProcessing=dataProcessing,
              processingOperation=processingOperation,
              scan=scan,
              scanOrigin=scanOrigin,
              precursorMz=precursorMz,
              maldi=maldi,
              peaks=peaks,
              sha1=sha1,
              .endDocument=calculateSha1Sums,
              .startElement=.startElement,
              .endElement=.endElement,
              .text=.text))
}
