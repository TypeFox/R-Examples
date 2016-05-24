#                   DF(+)/BDF(+) file header reading functions
#
# Purpose   :   Reading .edf(+) & .bdf file headers
#
# Copyright :   (C) 2015-2016, Vis Consultancy, the Netherlands
#               This program is free software: you can redistribute it and/or modify
#               it under the terms of the GNU General Public License as published by
#               the Free Software Foundation, either version 3 of the License, or
#               (at your option) any later version.
#
#               This program is distributed in the hope that it will be useful,
#               but WITHOUT ANY WARRANTY; without even the implied warranty of
#               MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#               GNU General Public License for more details.
#
#               You should have received a copy of the GNU General Public License
#               along with edfReader package for R.  If not, see <http://www.gnu.org/licenses/>.
#
# Acknowledgement :
#               This package has used code from:
#               - edf.R version 0.3 (27-11-2013) from Fabien Feschet,
#                 http://data-auvergne.fr/cloud/index.php/s/WYmFEDZylFWJzNs
#               - the work of Henelius Andreas as of July 2015, https://github.com/bwrc/edf
#
# History    :
#   Oct15 - Created
#   Feb16 - Version 1.0.0
#   Mar16 - Version 1.1.0 support for non-unique signals labels
#                         & sub second start data from data record
#   Apr16 - Version 1.1.1, no changes
# ------------------------------------------------------------------------------
#                                read EDF header
# ------------------------------------------------------------------------------
#' Read EDF(+) or BDF(+)  file header
#'
#' The function reads the header of an EDF (European Data Format) file, an EDF+ file,
#'    an BDF file, or an BDF+ file
#'
#' @param fileName The full path to the EDF(+)/BDF(+) file to be read.
#' @return an object of class ebdfHeader
#'
#' @section Details:
#' The object returned contains also an object of class ebdfSHeaders.
#' Both objects ebdfHeader and ebdfSHeaders have supporting S3 print and summary functions.\cr
#' For object details see the package vignette.
#'
#' @section Acknowledgement:
#'    This package has used code from:
#'    \itemize{
#'      \item edf.R version 0.3 (27-11-2013), http://feschet.fr/?p=11
#'      \item the work of Henelius Andreas as of July 2015, https://github.com/bwrc/edf
#'    }
#' @seealso
#'    \code{\link{edfReader}}, \code{\link{readEdfSignals}}\cr
#'    For the vignette use the console command:\cr
#'    \code{vignette('edfReaderVignette', package = "edfReader")}\cr
#'    or click on \code{Index} below.
#' @aliases readBdfHeader
#' @examples
#' # Examples from the vignette
#' libDir <- system.file ("extdata", package="edfReader")
#' # a continuous recording
#' CFile <- paste (libDir, '/edfPlusC.edf', sep='')
#' CHdr  <- readEdfHeader (CFile)
#' CHdr                   # print the header
#' summary (CHdr)         # print a header summary
#' str (CHdr)             # look to the details
#' CHdr$sHeaders          # print the signal headers
#' summary(CHdr$sHeaders) # print a signal headers summary
#' str(CHdr$sHeaders)     # look to the signal header details
#' # for a discontinuous recording
#' DFile <- paste (libDir, '/bdfPlusD.bdf', sep='')
#' # and proceed as above to read the header and to show the results
#' @export
readEdfHeader <- function (fileName) {
    # open the file
    if (file.exists (fileName)) inFile <- file(fileName, "rb", encoding="UTF-16LE")
    else stop (paste ("File '", fileName, "' doesn't exists.", sep=''))
    hdr <- list()
    class (hdr) <- "ebdfHeader"
    hdr$fileName    <- fileName
    hdr$fileType = ''
    version     <- readBin (inFile, integer(), n=1, size=1, signed=FALSE , endian="little")
    eVersion <- version
    if (version[1] == 48) {                                                     # ascii 48 = 0
        hdr$version <- '0'                                                      # always 0, see specs
        spaces <- readChar(inFile, 7, TRUE)                                     # read the other 7 bytes, should be 7 ' '
        if (spaces != '       ') stop (paste ("File '", fileName, "' is not an EDF or BDF file.", sep=''))
        hdr$fileType = 'EDF'
    } else if  (version[1] == 255) {                                            # possibly a bdf file
        hdr$fileType <- 'BDF'
        otherVersion <- gsub ("[[:space:]]*$", "", readChar(inFile, 7, TRUE))
        if (otherVersion == "BIOSEMI") {
            hdr$version <- paste('"', version,'"', otherVersion, sep='')
        } else {
            stop (paste ("File '", fileName, "' is not an EDF or BDF file.", sep=''))
        }
    } else stop (paste ("File '", fileName, "' is not an EDF or BDF file.", sep=''))

    hdr$patient     <- gsub ("[[:space:]]*$", "", readChar(inFile, 80, TRUE))   # local patient identification
    hdr$recordingId <- gsub ("[[:space:]]*$", "", readChar(inFile, 80, TRUE))   # local recording identification
    dateTimeSring   <- readChar(inFile, 16, TRUE)
    if (substr(dateTimeSring, 11, 11) == ':') {
        substr(dateTimeSring, 11, 11) <- '.'
        substr(dateTimeSring, 14, 14) <- '.'
    }
    hdr$startTime   <- strptime(dateTimeSring, format="%d.%m.%y %H.%M.%S")      # start date and time of recording
    # proposed?:       strptime(readChar(inFile, 16, TRUE), format="%d.%m.%y %H:%M:%S") # start date and time of recording
    hdr$headerLength<- as.integer(readChar(inFile, 8, TRUE))                    # number of bytes in header record
    hdr$reserved    <- gsub ("[[:space:]]*$", "", readChar(inFile, 44, TRUE))   # reserved
    hdr$nRecords    <- as.integer(readChar(inFile, 8, TRUE))                    # number of data records
    hdr$recordDuration <- as.numeric(readChar(inFile, 8, TRUE))                 # duration of a data record, in seconds
    hdr$nSignals    <- as.integer(readChar(inFile, 4, TRUE))                    # number of signals (ns) in data record

    # get the signal headers
    ns <- hdr$nSignals
    label       <- character (ns)
    for (i in 1:ns) label[i]            <- gsub ("[[:space:]]*$", "", readChar(inFile, nchars=16, useBytes=TRUE))
    transducerType <- character (ns)
    for (i in 1:ns) transducerType[i]   <- gsub ("[[:space:]]*$", "", readChar(inFile, nchars=80, useBytes=TRUE))
    physicalDim <- character (ns)
    for (i in 1:ns) physicalDim[i]      <- gsub ("[[:space:]]*$", "", readChar(inFile, nchars=8, useBytes=TRUE))
    physicalMin <- double (ns)
    for (i in 1:ns) physicalMin[i]      <-  as.numeric (readChar (inFile, nchars=8, useBytes=TRUE))
    physicalMax <- double (ns)
    for (i in 1:ns) physicalMax[i]      <- as.numeric (readChar (inFile, nchars=8, useBytes=TRUE))
    digitalMin <- integer (ns)
    for (i in 1:ns) digitalMin[i]        <- as.integer (readChar (inFile, nchars=8, useBytes=TRUE))
    digitalMax <- integer (ns)
    for (i in 1:ns) digitalMax[i]        <- as.integer (readChar (inFile, nchars=8, useBytes=TRUE))
    preFilter <- character (ns)
    for (i in 1:ns) preFilter[i]        <- gsub ("[[:space:]]*$", "", readChar(inFile, nchars=80, useBytes=TRUE))
    samplesPerRecord <- integer (ns)
    for (i in 1:ns) samplesPerRecord[i] <- as.integer (readChar (inFile, nchars=8, useBytes=TRUE))
    reserved <- character (ns)
    for (i in 1:ns) reserved[i]        <- gsub ("[[:space:]]*$", "", readChar(inFile, nchars=32, useBytes=TRUE))

    close(inFile)

    # derived attributes
    if (hdr$recordDuration)    sRate   <- samplesPerRecord / hdr$recordDuration
    else                       sRate   <- NA                                    # in EDF+ hdr$recordD may be zero
    sLength             <-   samplesPerRecord * hdr$nRecords
    hdr$recordedPeriod  <- hdr$recordDuration * hdr$nRecords

    hdr$sampleBits  <- ifelse (hdr$fileType == 'EDF', 16, 24)

    hdr$isPlus      <- FALSE                                                    # not an EDF+ file
    hdr$isContinuous<- TRUE                                                     # continuous recording, the default
    if (hdr$reserved == 'EDF+C' | hdr$reserved == 'BDF+C'| hdr$reserved == 'EDF+D'| hdr$reserved == 'BDF+D') {
        hdr$isPlus <- TRUE
        if (hdr$reserved == 'EDF+D' | hdr$reserved == 'BDF+D') hdr$isContinuous  <- FALSE
    }

    isAnnotation <- rep (FALSE, ns)
    if (hdr$isPlus) for (i in 1:ns) {
        if (label[i] == 'EDF Annotations' | label[i] == 'BDF Annotations') isAnnotation[i] <- TRUE
    }
    # physicaL =  gain * digital + offset  with:
    #     gain   = (physicalMax - physicalMin) / (digitalMax  - digitalMin)
    #     offset =  physicalMax - gain * digitalMax
    gain   <- rep (1, ns)
    offset <- rep (0, ns)
    for (i in 1:ns) if (!isAnnotation[i]) {
        if ((digitalMax[i]  != digitalMin[i]) & (physicalMax[i] != physicalMin[i])) {
            gain[i]   <- (physicalMax[i] - physicalMin[i]) / (digitalMax[i]  - digitalMin[i])
            offset[i] <- physicalMax[i] - gain[i] * digitalMax[i]
        }
    }
    # create signal data frame
    hdr$sHeaders <- data.frame(label=label, transducerType=transducerType, physicalDim=physicalDim,
                              physicalMin=physicalMin, physicalMax=physicalMax,
                              digitalMin=digitalMin, digitalMax=digitalMax,
                              preFilter=preFilter, samplesPerRecord=samplesPerRecord,
                              reserved=reserved, gain=gain, offset=offset, sRate=sRate,
                              isAnnotation=isAnnotation, sLength=sLength, stringsAsFactors=FALSE)
    signalNames <- label
    for (sn in 1:length(signalNames)) {
        if (signalNames[sn]=="") signalNames[sn] <- 'NA'                        # avoid empty names
    }
    for (sn in 1:length(signalNames)) {
        equals <- signalNames == signalNames[sn]
        if (sum(equals) > 1) {                                                  # distinguish duplicates
            suffix <- 1
            for (en in which(equals)) {
                # NOTE: This is not bullet proof "a-1, a, a" will result in "a-1, a-1, a-2".
                newName <- paste(signalNames[en], '-', suffix, sep='')
                while (sum(signalNames == newName)) {
                    suffix <- suffix +1
                    newName <- paste(signalNames[en], '-', suffix, sep='')
                }
                signalNames[en] <- paste(signalNames[en], '-', suffix, sep='')
                suffix <- suffix +1
            }
        }
    }

    #  if + file get secondfraction start form first data record
    hdr$startSecondFraction <- 0
    if (hdr$isPlus) {
        hdr$startSecondFraction <- getStartSecondFraction (hdr)
        hdr$startTime <- hdr$startTime + hdr$startSecondFraction
    }

    row.names (hdr$sHeaders) <- signalNames
    class (hdr$sHeaders) <- c("ebdfSHeaders", "data.frame")
    # perform a consistency check
    if (hdr$isPlus & !sum(hdr$sHeaders$isAnnotation)) cat ('ERROR: The annotation signal is missing.')
    hdr
}

getStartSecondFraction <- function (hdr) {
    StartSecondFraction <- 0
    nAnnots <- sum (hdr$sHeaders$isAnnotation)
    if (hdr$isPlus & nAnnots) {
        annotations1 <- which.max (hdr$sHeaders$isAnnotation)
        inFile <- file(hdr$fileName, "rb", encoding="UTF-16LE")
        readBin (inFile, 'raw', n=hdr$headerLength , size=1)
        samples             <- readNextDataRecord (hdr, inFile)
        StartSecondFraction <- getAnnoRecordStartRT (samples[[annotations1]])   # the subsecond start time is in the first TAL in the first annotation signal
        close(inFile)
    }
    StartSecondFraction
}
