##
## Copyright (c) 2010-2015, Brandon Whitcher
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * Neither the name of Rigorous Analytics Ltd. nor the names of
##       its contributors may be used to endorse or promote products
##       derived from this software without specific prior written
##       permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
## $Id: $
##

#' Read Single DICOM File
#'
#' All information, both header and image, is read into a list structure from a
#' DICOM file.
#'
#' A \code{while} loop is used to traverse the unknown number of DICOM header
#' fields contained in a single file.  Information contained in
#' \dQuote{sequences} may be included/excluded according to the logical
#' variable \code{skipSequence} (default = \code{TRUE}).
#'
#' A resursive implementation of the code breaks the DICOM file into segments
#' and calls itself to parse each segment.
#'
#' Strict adherence to the DICOM standard is not required.  Specifically,
#' content is allowed to start at the first byte and the four characters
#' \sQuote{DICM} are not required at bytes 129-132.
#'
#' @aliases parseDICOMHeader readDICOMFile
#' @param fname is the file name of the DICOM image (with suffix).
#' @param boffset is the number of bytes to skip at the beginning of the DICOM
#' file (default = \code{NULL} which lets the code determine the starting
#' point).
#' @param endian is the endian-ness of the file (default is \code{"little"}).
#' @param flipud is a logical variable for vertical flipping of the image
#' (default is \code{TRUE}).
#' @param skipSequence is a logical variable to skip all content contained in
#' SequenceItem tags (default = \code{TRUE}).
#' @param pixelData is a logical variable (default = \code{TRUE}) on whether or
#' not the PixelData should be read from the DICOM files.  This is useful when
#' one wants to gather the DICOM header information without loading the images.
#' @param warn is a number to regulate the display of warnings (default = -1).
#' See \code{options} for more details.
#' @param debug is a logical variable (default = \code{FALSE}) that regulates
#' to display of intermediate processing steps.
#' @param verbose is a logical variable (default = \code{FALSE}) that regulates
#' to display of intermediate processing steps.
#' @param rawString is a vector of \code{raw} values taken directly from the
#' DICOM file.
#' @param sq.txt is an character string (default = \dQuote{}) that indicates if
#' the DICOM header field is embedded within a sequence.
#' @return A list containing two elements: \describe{ \item{hdr}{all DICOM
#' header fields (with or without \dQuote{sequence} information).}
#' \item{img}{the \sQuote{image} information.} }
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{readDICOM}}
#' @references Whitcher, B., V. J. Schmid and A. Thornton (2011).  Working with
#' the DICOM and NIfTI Data Standards in R, \emph{Journal of Statistical
#' Software}, \bold{44} (6), 1--28.  \url{http://www.jstatsoft.org/v44/i06}
#'
#' Digital Imaging and Communications in Medicine (DICOM)\cr
#' \url{http://medical.nema.org}\cr
#' \url{http://en.wikipedia.org/wiki/Digital_Imaging_and_Communications_in_Medicine}
#' @keywords file
#' @examples
#'
#' x <- readDICOMFile(system.file("dcm/Abdo.dcm", package="oro.dicom"))
#' graphics::image(t(x$img), col=grey(0:64/64), axes=FALSE, xlab="", ylab="",
#'                 main="Abdo.dcm")
#'
#' x <- readDICOMFile(system.file("dcm/Spine1.dcm", package="oro.dicom"))
#' graphics::image(t(x$img), col=grey(0:64/64), axes=FALSE, xlab="", ylab="",
#'                 main="Spine1.dcm")
#'
#' @export readDICOMFile
readDICOMFile <- function(fname, boffset=NULL, endian="little", flipud=TRUE,
                          skipSequence=FALSE, pixelData=TRUE,
                          warn=-1, debug=FALSE) {
  ## Warnings?
  oldwarn <- getOption("warn")
  options(warn = warn)
  ##
  fsize <- file.info(fname)$size
  fraw <- readBin(fname, "raw", n=as.integer(fsize), endian=endian)
  if (is.null(boffset) && any(as.integer(fraw[1:128]) != 0)) {
    stop("Non-zero bytes are present in the first 128, please use\nboffset to skip the necessary number of bytes.")
  }
  skip128 <- fraw[1:128]
  skipFirst128 <- ifelse(any(as.logical(skip128)), FALSE, TRUE)
  if (debug) {
    cat("# First 128 bytes of DICOM header =", fill=TRUE)
    print(skip128)
  }
  DICM <- .rawToCharWithEmbeddedNuls(fraw[129:132]) == "DICM"
  if (debug) {
    cat("# DICM =", DICM, fill=TRUE)
  }
  #if (DICM) {
  #  if (.rawToCharWithEmbeddedNuls(fraw[129:132]) != "DICM") {
  #    stop("DICM != DICM")
  #  }
  #}
  dicomHeader <- sequence <- NULL
  seq.txt <- ""
  if (is.null(boffset)) {
    bstart <- 1 + ifelse(skipFirst128, 128, 0) + ifelse(DICM, 4, 0) # number of bytes to start
  } else {
    bstart <- boffset + 1
  }
  ## Call parseDICOMHeader() to recursively parse the DICOM header
  dcm <- parseDICOMHeader(fraw[bstart:fsize], seq.txt, endian=endian, verbose=debug)
  hdr <- as.data.frame(dcm$header, stringsAsFactors=FALSE)
  row.names(hdr) <- NULL
  names(hdr) <- c("group", "element", "name", "code", "length", "value", "sequence")
  ##
  if (dcm$pixel.data && pixelData) {
    if (debug) {
      cat("##### Reading PixelData (7FE0,0010) #####", fill=TRUE)
    }
    img <- parsePixelData(fraw[(bstart + dcm$data.seek):fsize], hdr, endian, flipud)
  } else {
    if (dcm$spectroscopy.data && pixelData) {
      if (debug) {
        cat("##### Reading SpectroscopyData (5600,0020) #####", fill=TRUE)
      }
      img <- parseSpectroscopyData(fraw[(bstart + dcm$data.seek):fsize], hdr, endian)
    } else {
      img <- NULL
    }
  }
  ## Warnings?
  options(warn = oldwarn)
  ##
  list(hdr = hdr, img = img)
}

.rawToCharWithEmbeddedNuls <- function(str.raw, to="UTF-8") {
  iconv(rawToChar(str.raw[str.raw != as.raw(0)]), to=to)
}
#' @rdname readDICOMFile
#' @export parseDICOMHeader
parseDICOMHeader <- function(rawString, sq.txt="", endian="little",
                             verbose=FALSE) {
  ##
  ## "The default DICOM Transfer Syntax, which shall be supported by
  ## all AEs, uses Little Endian encoding and is specified in Annex
  ## A.1." (PS 3.5-2004, page 38)
  ##
  ## PS 3.5-2004, Sect 7.1.2: Data Element Structure with Explicit VR
  ## Explicit VRs store VR as text chars in 2 bytes.
  ## VRs of OB, OW, SQ, UN, UT have VR chars, then 0x0000, then 32 bit VL:
  ##
  ## +-----------------------------------------------------------+
  ## |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 | 10 | 11 |
  ## +----+----+----+----+----+----+----+----+----+----+----+----+
  ## |<Group-->|<Element>|<VR----->|<0x0000->|<Length----------->|<Value->
  ##
  ## Other Explicit VRs have VR chars, then 16 bit VL:
  ##
  ## +---------------------------------------+
  ## |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |
  ## +----+----+----+----+----+----+----+----+
  ## |<Group-->|<Element>|<VR----->|<Length->|<Value->
  ##
  ## Implicit VRs have no VR field, then 32 bit VL:
  ##
  ## +---------------------------------------+
  ## |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |
  ## +----+----+----+----+----+----+----+----+
  ## |<Group-->|<Element>|<Length----------->|<Value->
  ##
  rawToHex <- function(bytes) {
    toupper(paste(rev(bytes), collapse=""))
  }
  is.item <- function (group, element) {
    group == "FFFE" && element %in% c("E000","E00D","E0DD")
  }
  ## Data files that are necessary to proceed
  #data(dicom.dic, package="oro.dicom", envir=environment())
  #data(dicom.VR, package="oro.dicom", envir=environment())
  ##
  strseek <- dseek <- 0
  dicomHeader <- NULL
  pixelData <- spectroscopyData <- FALSE
  while (! (pixelData || spectroscopyData) && strseek < length(rawString)) {
    ## rm(group, element, dictionaryIndex, dic, rawValue, VR, vr, value, length)
    group <- rawToHex(rawString[strseek + 1:2])
    element <- rawToHex(rawString[strseek + 3:4])
    if (! any(dictionaryIndex <- group == oro.dicom::dicom.dic$group & element == oro.dicom::dicom.dic$element)) {
      ## Private tag = Unknown
      dic <- list(group=group, element=element, code="UN", offset=1, name="Unknown")
    } else {
      dic <- oro.dicom::dicom.dic[dictionaryIndex, ]
    }
    if (verbose) {
      cat("#", group, element, dic$name, dic$code, sep="\t")
    }
    b56 <- .rawToCharWithEmbeddedNuls(rawString[strseek + 5:6])
    b78 <- readBin(rawString[strseek + 7:8], "integer", size=2, endian=endian)
    b47 <- readBin(rawString[strseek + 5:8], "integer", size=4, endian=endian)
    strseek <- strseek + 8
    if (b56 %in% c("OB","OW","SQ","UN","UT") && ! is.item(group, element)) {
      ## Explicit VR
      length <- readBin(rawString[strseek + 1:4], "integer", size=4, endian=endian)
      strseek <- strseek + 4
      if (b56 != "SQ") {
        rawValue <- rawString[strseek + 1:length]
      }
      vr <- b56
    } else {
      if (b56 %in% c("AE","AS","AT","CS","DA","DS","DT","FL","FD","IS","LO",
                     "LT","OF","PN","SH","SL","SS","ST","TM","UI","UL","US")
          && ! is.item(group, element)) {
        ## Explicit VR
        length <- b78
        if (dic$code != "SQ") {
          rawValue <- rawString[strseek + 1:length]
        }
        vr <- b56
      } else {
        ## Implicit VR
        length <- b47
        vr <- NULL
        if (is.item(group, element)) {
          length <- 0
          rawValue <- raw(0)
        } else {
          if (dic$code != "SQ") {
            rawValue <- rawString[strseek + 1:length]
          } else {
            vr <- "SQ"
          }
        }
      }
    }
    if (! is.null(vr)) {
      VR <- oro.dicom::dicom.VR[oro.dicom::dicom.VR$code == vr, ]
    } else if (dic$code != "UN") {
      VR <- oro.dicom::dicom.VR[oro.dicom::dicom.VR$code == dic$code, ]
    } else {
      VR <- oro.dicom::dicom.VR[oro.dicom::dicom.VR$code == "UN", ]
    }
    if (verbose) {
      cat("", VR$code, length, sep="\t")
    }
    if (sq.txt == "" && group == "7FE0" && element == "0010") {
      ## PixelData
      value <- "PixelData"
      pixelData <- TRUE
      dseek <- strseek
    } else {
      if (sq.txt == "" && group == "5600" && element == "0020") {
        ## SpectroscopyData
        value <- "SpectroscopyData"
        spectroscopyData <- TRUE
        dseek <- strseek + 4 # HACK: not sure why I need to skip an extra four bytes
      } else {
        if (VR$code %in% c("UL","US")) { # (VR$code == "UL" || VR$code == "US") {
          value <- readBin(rawValue, "integer", n=length/VR$bytes,
                           size=VR$bytes, signed=FALSE, endian=endian)
        } else if (VR$code %in% c("SL","SS")) { # (VR$code == "SL" || VR$code == "SS") {
          value <- readBin(rawValue, "integer", n=length/VR$bytes,
                           size=VR$bytes, signed=TRUE, endian=endian)
        } else if (VR$code %in% c("FD","FL")) { # (VR$code == "FD" || VR$code == "FL") {
          value <- readBin(rawValue, "numeric", n=length/VR$bytes,
                           size=VR$bytes, signed=TRUE, endian=endian)
        } else if (VR$code %in% c("OB","OW")) { # (VR$code == "OB" || VR$code == "OW") {
          value <- .rawToCharWithEmbeddedNuls(rawValue)
        } else if (VR$code == "SQ") {
          value <- "Sequence"
        } else {
          if (length > 0) {
            tmpString <- .rawToCharWithEmbeddedNuls(rawValue)
            tmpString <- sub(" +$", "", tmpString)     # remove white space at end
            tmpString <- gsub("[\\/]", " ", tmpString) # remove "/"s
            tmpString <- gsub("[\\^]", " ", tmpString) # remove "^"s
          } else {
            tmpString <- ""
          }
          value <- tmpString
        }
      }
    }
    if (verbose) {
      cat("", value, sq.txt, sep="\t", fill=TRUE)
    }
    dicomHeaderRow <- c(group, element, dic$name, dic$code, length, value, sq.txt)
    dicomHeader <- rbind(dicomHeader, dicomHeaderRow)
    if (group == "FFFE" && element == "E0DD") { # SequenceDelimitationItem
      dseek <- strseek
      break
    }
    if (VR$code == "SQ") {
      groupElement <- paste("(", group, ",", element, ")", sep="")
      if (length > 0) {
        ## Pass length of bytes provided explicitly by the sequence tag
        dcm <- parseDICOMHeader(rawString[strseek + 1:length],
                                paste(sq.txt, groupElement),
                                verbose=verbose)
      } else {
        ## Pass remaining bytes and look for SequenceDelimitationItem tag
        dcm <- parseDICOMHeader(rawString[(strseek + 1):length(rawString)],
                                paste(sq.txt, groupElement),
                                verbose=verbose)
        length <- dcm$data.seek
      }
      dicomHeader <- rbind(dicomHeader, dcm$header)
    }
    strseek <- strseek + ifelse(length >= 0, length, 0)
  }
  list(header = dicomHeader, pixel.data = pixelData, data.seek = dseek,
       spectroscopy.data = spectroscopyData)
}

#' Parse DICOM Pixel or Spectroscopy Data
#'
#' These subroutines process the information contained after the DICOM header
#' and process this information into an image (2D or 3D) or complex-valued
#' vector.
#'
#' A \code{while} loop is used to traverse the unknown number of DICOM header
#' fields contained in a single file.  Information contained in
#' \dQuote{sequences} may be included/excluded according to the logical
#' variable \code{skipSequence} (default = \code{TRUE}).
#'
#' A resursive implementation of the code breaks the DICOM file into segments
#' and calls itself to parse each segment.
#'
#' @aliases parsePixelData parseSpectroscopyData
#' @param rawString is a vector of \code{raw} values taken directly from the
#' DICOM file.
#' @param hdr is the list object of DICOM header information.
#' @param endian is the endian-ness of the file (default is \code{"little"}).
#' @param flipupdown is a logical variable for vertical flipping of the image
#' (default is \code{TRUE}).
#' @return A list containing two elements: \describe{ \item{hdr}{all DICOM
#' header fields (with or without \dQuote{sequence} information).}
#' \item{img}{the \sQuote{image} information.} }
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{parseDICOMHeader}}, \code{\link{readDICOMFile}}.
#' @references Digital Imaging and Communications in Medicine (DICOM)\cr
#' \url{http://medical.nema.org}\cr
#' \url{http://en.wikipedia.org/wiki/Digital_Imaging_and_Communications_in_Medicine}
#' @source See references.
#' @keywords file
#' @export parsePixelData
parsePixelData <- function(rawString, hdr, endian="little", flipupdown=TRUE) {
  rows <- as.numeric(with(hdr, value[name == "Rows" & sequence == ""]))
  columns <- as.numeric(with(hdr, value[name == "Columns" & sequence == ""]))
  bytes <- as.numeric(with(hdr, value[name == "BitsAllocated" & sequence == ""])) / 8
  length <- as.numeric(with(hdr, length[name == "PixelData" & sequence == ""]))
  if (length <= 0) {
    guess <- 1
    stop(paste("Number of bytes in PixelData not specified; guess =", guess))
  }
  pixelRepresentation <- as.numeric(with(hdr, value[name == "PixelRepresentation" & sequence == ""]))
  signed <- ifelse(pixelRepresentation == 1, TRUE, FALSE)
  imageData <- readBin(rawString[1:length], "integer", n=length, size=bytes,
                       signed=signed, endian=endian)
  if (length == rows * columns * bytes) { # 2D PixelData
    imageData <-  t(matrix(imageData[1:(columns * rows)], columns, rows))
    if (flipupdown) {
      imageData <- imageData[rows:1, ]
    }
  } else { # 3D PixelData
    k <- length / rows / columns / bytes
    if (k == trunc(k)) {
      warning("3D DICOM file detected!")
      samplesPerPixel <- with(hdr, value[grepl("SamplesPerPixel", name, ignore.case=TRUE) & sequence == ""])
      if (samplesPerPixel == "1") { # && planarConfiguration == 0) {
        imageData <- array(imageData[1:(columns * rows * k)], c(columns, rows, k))
        imageData <- aperm(imageData, c(2,1,3))
        if (flipupdown) {
          imageData <- imageData[rows:1, , ]
        }
      } else {
        planarConfiguration <- with(hdr, value[grepl("PlanarConfiguration", name, ignore.case=TRUE) & sequence == ""])
        if (planarConfiguration == "0")
          stop("Color channels are interlaced")
      }
    } else {
      stop("Number of bytes in PixelData does not match dimensions of image")
    }
  }
  return(imageData)
}
#' @rdname parsePixelData
#' @export parseSpectroscopyData
parseSpectroscopyData <- function(rawString, hdr, endian="little") {
  numberOfFrames <- as.numeric(with(hdr, value[name == "NumberOfFrames" & sequence == ""]))
  rows <- as.numeric(with(hdr, value[name == "Rows" & sequence == ""]))
  columns <- as.numeric(with(hdr, value[name == "Columns" & sequence == ""]))
  dataPointRows <- as.numeric(with(hdr, value[name == "DataPointRows" & sequence == ""]))
  dataPointColumns <- as.numeric(with(hdr, value[name == "DataPointColumns" & sequence == ""]))
  dataRepresentation <- with(hdr, value[name == "DataRepresentation" & sequence == ""])
  nComponents <- ifelse(dataRepresentation == "COMPLEX", 2, 1)
  whichComponent <- ifelse(nComponents > 1 && dataRepresentation == "IMAGINARY", 2, 1)
  valuesPerFrame <- columns * rows * dataPointRows * dataPointColumns
  bytes <- as.numeric(with(hdr, value[name == "BitsAllocated" & sequence == ""])) / 8
  length <- valuesPerFrame * numberOfFrames
  imageData <- readBin(rawString[1:length], "numeric", n=length, size=4, endian=endian)
  odd <- seq(1, length, by=2)
  even <- seq(2, length, by=2)
  complex(real=imageData[odd], imaginary=imageData[even])
}
