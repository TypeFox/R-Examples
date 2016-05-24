#                      EDF(+)/BDF(+) file signals reader functions
#
# Purpose   :   Reading signals in  .edf(+)/.bdf(+) files
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
# History    :
#   Oct15 - Created
#   Feb16 - Version 1.0.0
#   Mar16 - Version 1.1.0 with:
#               - support for multiple annotation signals; multiple annatations per TAL
#               - annotations data.frame with only one annotation per row (resolves a bug)
#               - support sub second start specification in first data record (in + file)
#               - signal names copied from sHeaders
#  Apr16 - Version 1.1.1, no changes
#
# Suffixes (not always used):
# HRT       : header reference time, i.e. time in sec relative to start of recording as specified in the header.
#             one unaligned value for all signals
# RRT       : recording reference time, i.e. HRT + the subsecond specified as start for the first record, if present.
#             one unaligned value for all signals
#             the 'form' and 'till' parameters are interpreted as RRTs
# ART       : aligned recording timne, i.e. aligned to the signals sample rate started at RRT=0
#             The difference between RRT and ART is relevant only if the RRT for the start of a data record is not
#             a muliple of the record duration as specified in the header, i.e. for +D files only.
# RS        : Sample number relative signal sampling startin at RRT=0
#
#             (but note that that the actual start may differ if the first TAL in the first annotation record != 0)
# DS and DT : Data record Sample number en Time in seconds, i.e. relative to the start of a data record
# FS and FT : Fragment Sample number en Time in seconds, i.e. relative to the start of signal fragment
# L and D   : Length (number of samples) and duration (in seconds)
#
# prefixes indicate objects (not always used)
# wRec      : the whole recording
# dRec      : a data record  (using URT)
# recS      : a data record signal (using ART)
# sel       : the selection made with the from and till parameters
#
# ------------------------------------------------------------------------------
#                        read ordinary / annotation signals
# ------------------------------------------------------------------------------
#' Reads signals from an EDF(+)/BDF(+) file
#'
#' The function reads ordinary or annotation signals from an EDF(+)/BDF(+) file.
#'
#' The list of signals returned is of class ebdfSignals and a single signal object is of one of the
#' following classes:
#' \itemize{
#'   \item ebdfASignal, for an annotation signal
#'   \item ebdfFSignal, for a fragmented ordinary signal
#'   \item ebdfCSignal, for a continuous ordinary signal (possible supplemented with NA values)
#' }
#' All classes have supporting print and summary functions.
#' For object details see the package vignette.
#' @param hdr An ebdfHeader object read with the readEdfHeader() function.
#' @param signals a vector with one or more of the following signal designations:
#'  'All' (default), to include all signals;
#'  'Ordinary', to include all ordinary signals;
#'  'Annotations', to include all annotation signals;
#'  signal labels and/or signal names
#'  signal numbers (numeric or as character).
#' @param from numeric, the time in seconds from which the signals shall be read.
#' @param till numeric, the time in seconds till which the signals shall be read.
#'   The value may exceed the total duration of the recoding.
#' @param physical logical, if TRUE (the default) digital samples values are mapped to their physical values,
#'   If not, the digital values are returned.
#' @param fragments logical, if TRUE discontinuously recorded signals are stored as a list of continuous
#'   fragments. If FALSE (the default), a signal is stored as one numeric vector with NA values filling the gaps.
#' @param recordStarts logical, if TRUE the empty annotations with the data record start time will be
#'   included in the list of annotations.
#'   If FALSE (the default), they will be omitted.
#' @param mergeASignals logical, if TRUE all annotations will wil merged into one data frame. If FALSE
#'   there will be one data frame per annotation signal.
#' @param simplify logical, if TRUE and if there is only one signal read, the signal itsels is returned
#'   (in stead of a list with that signal as the only one object).
#'   If FALSE, this simplification is not performed.
#' @return Either a list of one or more signals or a single signal.
#' @section Note:
#'   For ordinary signals the from and till parameters are interpreted as [from, till).
#'   For for details see the package vignette.
#' @section Acknowledgement:
#'    This package has used code from:
#'    \itemize{
#'      \item edf.R version 0.3 (27-11-2013), http://feschet.fr/?p=11
#'      \item the work of Henelius Andreas as of July 2015, https://github.com/bwrc/edf
#'    }
#' @seealso
#'    \code{\link{edfReader}}, \code{\link{readEdfHeader}}\cr
#'    For the vignette use the console command:\cr
#'    \code{vignette('edfReaderVignette', package = "edfReader")}\cr
#'    or click on \code{Index} below.
#' @aliases readBdfSignals
#' @examples
#' # Examples from the vignette
#' libDir <- system.file ("extdata", package="edfReader")
#' # a continuous recording
#' CFile <- paste (libDir, '/edfPlusC.edf', sep='')
#' CHdr  <- readEdfHeader (CFile)
#' CSignals <- readEdfSignals (CHdr)    # to read all signals
#' # read 3 differently designated signals from 5.1 till 18 seconds period
#' someCSignalsPeriod <- readEdfSignals (CHdr, signals=c(3, "5", "sine 8.5 Hz"), from=5.1, till=18)
#' someCSignalsPeriod                           # print the signals
#' summary(someCSignalsPeriod)                  # print singals summary
#' someCSignalsPeriod$`sine 8.5 Hz`             # print the `sine 8.5 Hz` signal
#' summary(someCSignalsPeriod$`sine 8.5 Hz`)    # print a `sine 8.5 Hz` signal summary
#' str(CSignals$`sine 8.5 Hz`)                  # look to the details
#' # a discontinuous recording
#' DFile <- paste (libDir, '/edfPlusD.edf', sep='')
#' DHdr  <- readEdfHeader (DFile)
#' DSignals <- readEdfSignals (DHdr, fragments=TRUE)    # to read all signals
#' DSignals$`sine 8.5 Hz`                       # print fragmented signal
#' summary (DSignals$`sine 8.5 Hz`)             # print fragmented signal summary
#' str(DSignals$`sine 8.5 Hz`)                  # look to the details
#' @export
readEdfSignals <- function (hdr, signals='All', from=0, till=Inf, physical=TRUE,
                            fragments=FALSE, recordStarts=FALSE, mergeASignals=TRUE,
                            simplify=TRUE) {

    # check signals and get the indices
    idxAndErr <- edfProcessSignalDesignations (hdr, signals)
    if (is.null(idxAndErr)) stop ("Illegal 'signals' parameter")
    errs <- sum (idxAndErr$errors)
    if (errs) {
        if (errs==1)    pre <- 'Unkown signal designation'
        else            pre <- 'Unkown signal designations'
        stop (paste (pre, signals[idxAndErr$errors]))
    }
    if (length(idxAndErr$outOfBound)) {
        if (length(idxAndErr$outOfBound) == 1)  pre <- "Signal number out of bound:"
        else                                    pre <- "Signal numbers out of bound:"
        stop (paste(pre, paste(idxAndErr$outOfBound, collapse = " ")))
    }

    if (physical) {                                                             # conversion mut be possible
        digitalOk    <- hdr$sHeaders$digitalMin  <  hdr$sHeaders$digitalMax
        physicalOk   <- hdr$sHeaders$physicalMin != hdr$sHeaders$physicalMax
        digitalErr   <- sum (!digitalOk)
        physicalErr  <- sum (!physicalOk)
        if (digitalErr | physicalErr) {
            if (digitalErr & physicalErr) {
                msg <- "Illegal digital/physical min/max, use physical=FALSE"
            } else if (digitalErr) {
                msg <- "Illegal digital min/max, use physical=FALSE"
            } else if (physicalErr) {
                msg <- "Illegal physical min/max, use physical=FALSE"
            }
            stop (msg)
        }
    }

    # check for annotation and non continuous signals in EDF+ and BDF+
    isPlus       <- hdr$isPlus
    isContinuous <- hdr$isContinuous
    sgn          <- idxAndErr$signals                                           # signals to be returned
    oSignals     <- sgn & !hdr$sHeaders$isAnnotation                            # the OSignals
    annotations1 <- which.max (hdr$sHeaders$isAnnotation)
    aSignals     <- sgn & hdr$sHeaders$isAnnotation
    aSignals[annotations1] <- TRUE                                              # contains record start times

    if (isPlus & !sum(hdr$sHeaders$isAnnotation)) {
        fp <- paste (hdr$fileType, '+', sep='')
        cat (fp, 'ERROR: required annotation signal is missing, signals will be retrieved as for an',
             hdr$fileType, 'file.\n')
        isPLus <- FALSE
        isContinuous <- TRUE
    }

    # open file; if +C get start time first record and reopen; skip the header  NOTE : seek is not advised for Windows
    if (file.exists (hdr$fileName)) inFile <- file(hdr$fileName, "rb", encoding="UTF-16LE")
    else stop (paste ("File '", hdr$fileName, "' doesn't exists.", sep=''))
    readBin (inFile, 'raw', n=hdr$headerLength , size=1)                        # the header is not needed;

    # check and justify 'from' and 'till'
    if (from > till) {
        stop ("Illegal from-till range.")
    }
    wRecFromHRT <- hdr$startSecondFraction                                      # as.numeric (hdr$startTime - trunc (hdr$startTime, 'sec'))
    if (isContinuous & from > hdr$recordedPeriod) {
        stop ("From after the continuously recorded period.")
    }

    # initialise signals
    signals <- vector (mode='list', length=hdr$nSignals)                        # new list with a element per signal
    names (signals) <- row.names (hdr$sHeaders)                                 # signalNames
    nextInCSignal   <- integer (length=hdr$nSignals) + 1                        # used for ordinary continuous signals only

    # intialise signal boundary parameters
    maxTErr         <- 5 * .Machine$double.eps                                  # normally 5 * 2.220446e-16; to avoid rounding errors when ceiling
    # aligned start of recording per signal
    wRecFromART     <- 0                                                        # aligned first recorded sample

    # aligned selection per signal
    selFromRRT    <- max (0, from)
    selTillRRT    <- till
    if (isContinuous) {
        selTillRRT    <- min (hdr$recordDuration * hdr$nRecords, till)
    }
    selFromRS     <- ceiling (hdr$sHeaders$sRate * (selFromRRT - maxTErr)) +1   # aligned first sample in selection per signal
    selFromART    <- (selFromRS-1) / hdr$sHeaders$sRate                         # aligned start of selection per signal
    selTillRS     <- ceiling (hdr$sHeaders$sRate * (selTillRRT - maxTErr))      # aligned last sample in selection per signal
    selTillART    <- (selTillRS-1) / hdr$sHeaders$sRate                         # aligned end of selection per signal

    # initialise data objects to be read
    for (sn in 1:hdr$nSignals) if (sgn[sn]) {                                   # => NULL if not in sgn
        signals[[sn]]$startTime     <- hdr$startTime
        signals[[sn]]$signalNumber  <- sn
        signals[[sn]]$label         <- hdr$sHeaders$label[sn]
        signals[[sn]]$isContinuous  <- isContinuous
        signals[[sn]]$isAnnotation  <- hdr$sHeaders$isAnnotation[sn]
        if (oSignals[sn]) {
            signals[[sn]]$recordedPeriod<- hdr$recordedPeriod
        }
        signals[[sn]]$totalPeriod   <- as.numeric(NA)
        signals[[sn]]$from          <- from
        signals[[sn]]$till          <- till
        if (oSignals[sn]) {
            range <- paste (hdr$sHeaders$physicalMin[sn], " : ", hdr$sHeaders$physicalMax[sn],
                            ' ', hdr$sHeaders$physicalDim[sn], sep = '')
            signals[[sn]]$start         <- selFromART[sn]
            signals[[sn]]$fromSample    <- selFromRS[sn]
            signals[[sn]]$transducerType<- hdr$sHeaders$transducerType[sn]
            signals[[sn]]$sampleBits    <- hdr$sampleBits
            signals[[sn]]$sRate         <- hdr$sHeaders$sRate[sn]
            signals[[sn]]$range         <- range
            signals[[sn]]$preFilter     <- hdr$sHeaders$preFilter[sn]
        }
        if (aSignals[sn]) {
            signals[[sn]]$annotations <- vector("list", length=hdr$nRecords)    # initialise for raw data per record
            class(signals[[sn]]) <- 'ebdfASignal'
        } else if (isContinuous) {
            for (sn in 1:hdr$nSignals) if (oSignals[sn]) {
                l       <- max (0, selTillRS[sn] - selFromRS[sn] + 1)
                signals[[sn]]$signal <- integer (length=l)                      # initialise with 0 signal
                if (l==0) {
                    signals[[sn]]$start         <- as.numeric (NA)
                    signals[[sn]]$fromSample    <- as.numeric (NA)
                }
                class(signals[[sn]]) <- 'ebdfCSignal'
            }
        } else {   # +D signal
            signals[[sn]]$rFragments <- vector("list", length=hdr$nRecords)     # initialise for raw data per record
            class(signals[[sn]]) <- 'ebdfFSignal'
        }
    }

    # read
    # cat ("Reading file", hdr$fileName, "signals, please wait.\n")
    for (rn in 1:hdr$nRecords) {
        samples <- readNextDataRecord (hdr, inFile)                             # read next data record
        # get the start data
        dRecFromRRT <- hdr$recordDuration * (rn-1)                              # if continuous
        dRecFromRS  <- hdr$sHeaders$samplesPerRecord * (rn-1) + 1               # idem
        dRecFromART <- rep (dRecFromRRT, times= hdr$nSignals)                   # idem
        if (isPlus) {
            dRecFromHRT <- getAnnoRecordStartRT (samples[[annotations1]])       # stated record start (in annotations1)
            if (hdr$isContinuous) {
                maxStartTDiff   <- 1e-10                                        # more tolerant, only used for a warning
                if ( abs(dRecFromHRT - dRecFromRRT -  wRecFromHRT) > maxStartTDiff) {
                    cat ("WARNING: Ambiguous data record start:\n",
                         "- start according to annotation signal:", dRecFromHRT - wRecFromHRT, '\n',
                         "- start based on recordDuration:", dRecFromRRT, '\n')
                }
            } else {                                                            # if not conintuous
                dRecFromRRT <- dRecFromHRT - wRecFromHRT
                dRecFromRS  <- ceiling (hdr$sHeader$sRate * (dRecFromRRT - maxTErr)) + 1
                dRecFromART <- (dRecFromRS-1) / hdr$sHeaders$sRate
            }
        }

        # (pre)process the signals to return
        for (sn in 1:hdr$nSignals) if (sgn[sn]) {
            if (hdr$sHeaders$isAnnotation[sn]) {
                signals[[sn]]$annotations[[rn]] <- samples[[sn]]
            } else {                                                            # an ordinary signal
                # set read boundaries for this record
                sLength <- hdr$sHeader$samplesPerRecord[sn]                     # signal samples per record
                skipRecord      <- FALSE
                wholeRecord     <- TRUE
                drsFromRS       <- dRecFromRS[sn]                               # dataRecordSignal
                drsFrommART     <- dRecFromART[sn]
                if (selFromRS[sn] > 0 | till < Inf) {
                    fromDS      <- max (1, selFromRS[sn] -  drsFromRS + 1)
                    tillDS      <- min (sLength, selTillRS[sn] - drsFromRS + 1)
                    skipRecord  <- (fromDS > sLength) | (tillDS <= 0)
                    wholeRecord <- (fromDS == 1) & (tillDS == sLength)
                }

                # convert sample values
                if (!skipRecord) {
                    # FIRST: for bdf, convert to integer; then samples[[sn]] will get the right length too
                    if (hdr$fileType == 'BDF') samples[[sn]] <- int1sToInt3s (samples[[sn]])
                    # truncate, if required
                    if (!wholeRecord) {
                        samples[[sn]] <- samples[[sn]] [fromDS:tillDS]
                    }
                    # convert digital to physical, if required
                    if (physical) {
                        samples[[sn]] <- hdr$sHeaders$gain[sn] * samples[[sn]] + hdr$sHeaders$offset[sn]
                    }
                    # add the samples to signals / rFragments
                    if (isContinuous) {
                        lastOne <- nextInCSignal[sn] + length (samples[[sn]]) -1
                        signals[[sn]]$signal[nextInCSignal[sn]:lastOne] <- samples[[sn]]  # copy to signals[[sn]]$signal
                        nextInCSignal[sn] <- lastOne + 1
                    } else {                                                    # not a continuous signal
                        fsFromRS    <- drsFromRS + fromDS - 1
                        fsFromART   <- (fsFromRS-1) / hdr$sHeaders$sRate[sn]
                        fsFromRRT   <- max (dRecFromRRT, selFromRRT)
                        signals[[sn]]$rFragments[[rn]] <- list (fsFromART=fsFromART, fsFromRS=fsFromRS,
                                                                fsFromRRT=fsFromRRT,
                                                                drsFromRS=drsFromRS,
                                                                signal=samples[[sn]])
                    }
                }
            }
        }
    }
    close(inFile)

    wRecTillRRT <- dRecFromRRT + hdr$recordDuration
    wRecTillRS  <- dRecFromRS + hdr$sHeader$samplesPerRecord - 1                # recordStartRS = aligned first sample last record
    wRecTillART <- wRecTillRS / hdr$sHeaders$sRate                              # = total recording period
    # save total period, recorded period (which may be part of the total period for EDF-D files)
    for (sn in 1:hdr$nSignals) if (sgn[sn]) {
        signals[[sn]]$totalPeriod       <- wRecTillRRT                          # including gaps
#        signals[[sn]]$recordedPeriod    <- wRecTillART
    }

    # process annotations
    for (sn in 1:hdr$nSignals) if (sgn[sn] & hdr$sHeaders$isAnnotation[sn]) {
        signals[[sn]]$annotations <- edfProcesAnnotations (
            hdr=hdr, ASignal=signals[[sn]],
            isFirstASignal = sn==annotations1, recordStarts=recordStarts
        )
    }
    # merge ASignals if requested
    if (mergeASignals) {
        annoSL <- hdr$sHeaders$isAnnotation & sgn
        if (sum(annoSL) > 1) {
            annoSN  <- which (annoSL)
            annoSN1 <- which.max (annoSL)
            signals[[annoSN1]]  <- doMergeASignals (signals[annoSN], annoSN)
            annoSL[annoSN1]     <- FALSE
            sgn     <- sgn & !annoSL

            names           <- names (signals)
            names[annoSN1]  <- substr (names[annoSN1], 1, 15)                   # remove suffix
            names(signals)  <- names
        }
    }
    # sort annotations
    annoSL <- hdr$sHeaders$isAnnotation & sgn
    if (sum(annoSL)) {
        for (sn in 1:hdr$nSignals) if (annoSL[sn]) {
            annots  <- signals[[sn]]$annotations
            annots  <- annots[order(annots$onset),]
            signals[[sn]]$annotations <- annots
        }
    }

    # process D signals
    if (!isContinuous & sum(oSignals)) {                                        # i.e.  +D file with ordinary siganals
        for (sn in 1:hdr$nSignals) if (oSignals [sn]) {

            # remove empty signal rFragments
            usedRFragments <- !sapply(signals[[sn]]$rFragments, is.null)
            signals[[sn]]$rFragments <- signals[[sn]]$rFragments[usedRFragments]

            # and check that no empty signals remain
            signalsL <- sapply (signals[[sn]]$rFragments, function(F) {length (F$signal)})
            noSignal <- signalsL == 0
            if (sum(noSignal)) cat ("readEdfSignals ERROR: Fragments with no signal !!!!!!!!!\n")

            if (!fragments) { # create a continuous signal with NAs
                till <- min (wRecTillRS[sn], selTillRS[sn])
                signals[[sn]] <- fragmentsToSignal (signals[[sn]], selFromRRT, selFromRS[sn], till)
            } else {  # fragmented, concatenate contiguous DSignal rFragments
                signals[[sn]]$fragments <- concatenateFragments (signals[[sn]]$rFragments, hdr$sHeaders$samplesPerRecord[sn])
            }
            signals[[sn]]$rFragments <- NULL
        }
    }

    # return the signals requested
    signals <- signals [sgn]
    class (signals) <- c("ebdfSignals", "list")
    if (simplify & length (signals)==1) signals <- signals[[1]]
    signals
}

# ------------------------------------------------------------------------------
#                        supplementary functions
# ------------------------------------------------------------------------------
# Maps signal designations to signal indices
#
# The function maps a list a signal designations into an signal index vector.
# It also reports illegal designations.
#
# @param edf The edf object with the EDF file header
# @param signals='All' which may have one of the following values:
#  - 'All', to include all signals
#  - 'Signals', to include all signal signals only
#  - 'Annotations', to include all annotation signals only
#  - A vector with signal designations consisting of labels and/or signal numbers (numeric or as character)
# @return A list with two logical vectors:
#   - One to indicate the selected signals
#   - One to indicate erroneous signal values
#
# @keywords internal
edfProcessSignalDesignations <- function (hdr, signals='All') {
    errorsL     <- rep (FALSE, length(signals))
    indicesL    <- rep (FALSE, hdr$nSignals)
    outOfBound  <- integer (length=0)
    for (sn in 1:length(signals)) {
        if      (signals[sn] == 'All')           indicesL[] <- TRUE
        else if (signals[sn] == 'Ordinary')      indicesL <- indicesL | !hdr$sHeaders$isAnnotation
        else if (signals[sn] == 'Annotations')   indicesL <- indicesL |  hdr$sHeaders$isAnnotation
        else  {
            idx <- suppressWarnings (as.integer(signals[sn]))                   # try a (coerced) integer
            if (!is.na(idx))
                if (idx <= hdr$nSignals) indicesL [idx] <- TRUE
                else                     outOfBound     <- c(outOfBound, idx)
            else {
                labels  <- hdr$sHeaders$label
                labelsL <- labels == signals[sn]
                rNames  <- row.names(hdr$sHeaders)
                rNamesL <- rNames == signals[sn]
                someFound <- sum (labelsL) + sum (rNamesL)
                if (someFound) {
                    indicesL <- indicesL | labelsL |  rNamesL
                } else{
                    errorsL [sn] <- TRUE
                }
            }
        }
    }
    list (signals=indicesL, errors=errorsL, outOfBound=outOfBound)
}

readNextDataRecord <- function (hdr, inFile) {
    samples     <- vector (mode='list', length=hdr$nSignals)
    sampleSize  <- hdr$sampleBits / 8
    for (sn in 1:hdr$nSignals) {                                                # read all signals (i.e. don't use a seek)
        n  <- hdr$sHeaders$samplesPerRecord[sn]
        # read all record data
        if (sampleSize ==2 & !hdr$sHeaders$isAnnotation[sn]) {                  # ordinary 16 bits signal
            samples[[sn]] <- readBin (inFile, integer(), n=n, size=sampleSize,   signed=TRUE, endian="little")
        } else {                                                                # annotation or 24 bits
            samples[[sn]] <- readBin (inFile, integer(), n=n*sampleSize, size=1, signed=TRUE, endian="little")
        }
    }
    samples
}

int1sToInt3s <- function (int1s) {                                              # int1s: array of *signed* 1 byte integers
    ni <- length (int1s) / 3
    samples <- numeric(length = ni)
    for (i in 1:ni) {
        j <- i * 3                                                              # the last and most signifcant sample byte
        if (int1s[j-1] < 0)  int1s[j-1] <- int1s[j-1] + 256                     # correct the less significant int1s
        if (int1s[j-2] < 0)  int1s[j-2] <- int1s[j-2] + 256
        samples [i] <- int1s[j]*65536 + int1s[j-1]*256 + int1s[j-2]
    }
    samples
}
# ------------------------------------------------------------------------------
#                        Annotation signal functions
# ------------------------------------------------------------------------------
# get the start time for a data record
#
# The function retieves the first tal from (the first) annotation signal
#
# @param annotationSignal The first annoation signal from the data record
# @return The start time in seconds
# @keywords internal
getAnnoRecordStartRT <- function (annotationSignal) {
    endings <- which(annotationSignal==0)                                       # the TAL endings / separators 0
    tal1    <- annotationSignal [1:(endings[1]-1)]
    pt      <- parseTal (tal1)
    pt$onset
}

# Copies all annotations into a data frame
#
# The function copies all annotation in an 'EDF Annotations' signal
#     into a data frame.
#
# @param hdr The hdr object with the EDF / BDF file header
# @param signal The index of the signal with the annotations
# @return A list with the following values: record, onset, duration,
#      isRecordStart, annotations
#
# @keywords internal
edfProcesAnnotations <- function (hdr, ASignal, isFirstASignal, recordStarts) {
    nAnnots <- 0
    nTals   <- 0
    from <- ASignal$from - .Machine$double.eps                                  # to mitigate rounding errors
    till <- ASignal$till + .Machine$double.eps
    annots <- ASignal$annotations

    for (rn in 1:hdr$nRecords) {                                                # strip & count for each record
        rnAnnots <- annots[[rn]]
        # remove trailing zero's
        i <- length(rnAnnots)
        while (i>0 && rnAnnots[i]==0) i <- i-1  #  !!!! or check for last '20' !!!!!!!!!!
        rnAnnots  <- rnAnnots[1:(i+1)]                                          # keep last TAL delimiter
        # add the number of TALs endings, i.e. the number of '0's
        if (length(rnAnnots) > 1) {
            nTals   <- nTals   + sum (rnAnnots==0)
            nAnnots <- nAnnots + sum (rnAnnots==20)                             # still needs adjustment
        }
        annots[[rn]] <- rnAnnots
    }
    # adjust nAnnots
    nAnnots <- nAnnots - nTals                                                  # per TAL: annots in between '20's
    if (isFirstASignal & !recordStarts) nAnnots <- nAnnots -hdr$nRecords        # ommit the recordstart times

    # create empty data frame (possibly too  long because of a from - till range)
    annotations <- data.frame(record=integer(nAnnots), onset=numeric(nAnnots),
                              duration=numeric(nAnnots), isRecordStart=logical(nAnnots),
                              annotation=character(nAnnots), stringsAsFactors=FALSE
                              )
    nextAnnon <- 1
    for (rn in 1:hdr$nRecords) {                                                # for each record
        rnAnnots <- annots[[rn]]
        endings <- which (rnAnnots==0)                                          # the TAL endings / separators 0
        ending1 <- endings[1]
        fromChar<- 1
        if (length(rnAnnots) > 1) for (it in endings) {
            tal  <- rnAnnots[fromChar:(it-1)]                                   # tal which trailing 20
            fromChar <- it + 1
            pt   <- parseTal (tal)
            if (from <=pt$onset & pt$onset < till) {
                anns    <- pt$annotations
                for (ia in 1:length(anns)) {
                    isRecordStart <- isFirstASignal & it==ending1 & ia ==1
                    if (isRecordStart) {                                        # the firts annoatation must be empty
                        if (anns[ia]!="") {
                            cat ('Illegal annotation signal, the start time annotation must be empty\n')
                        }
                    }
                    if (recordStarts || !isRecordStart) {
                        annotations$record[nextAnnon]       <- rn
                        annotations$onset[nextAnnon]        <- pt$onset
                        annotations$duration[nextAnnon]     <- pt$duration
                        annotations$isRecordStart[nextAnnon]<- isRecordStart
                        annotations$annotation[nextAnnon]   <- anns[ia]
                        nextAnnon <- nextAnnon +1
                    }
                }
            }
        }
    }
    if (nextAnnon == 1) {
        # remove if record == 0, i.e. remove the empty row
        annotations <- annotations[annotations$record, ]
    } else if (nextAnnon - 1 < nAnnots) annotations <- annotations[1:(nextAnnon - 1),]
    annotations
}

#  parses a TAL
#
# The function parses a TAL (Time-stamped Annotations List)
#
# @param tal A raw TAL from the edf / bdf file
# @param isRecordStart True if the TAL is the first one in a record, FALSE if not
# @return A list with the following values: onset, duration, annotations
# @keywords internal
parseTal <- function (tal) {
    endings <- which(tal==20)                                                   # locate 'phrase' trailing delimiters '20'
    onsetPlusDuration <- tal[1:(endings[1]-1)]                                  # without trailing 20
    di <- which(onsetPlusDuration==21)                                          # locate start delimiter for 'duration'
    if (length(di)) {
        onset   <- as.numeric (intToUtf8 (onsetPlusDuration[1:(di-1)]))
        duration<- as.numeric (intToUtf8 (onsetPlusDuration[(di+1):length(onsetPlusDuration)]))
    } else {
        onset   <- as.numeric (intToUtf8 (onsetPlusDuration))
        duration<- NA
    }
    l <-length (endings)
    if (l-1)    annotations <- character(l-1)                                   # l-1 = number of annotations
    else        annotations <- ''
    if (l>1) for (i in 2:l) {
        if (endings[i-1] == endings[i]-1) annotations[i-1] <- ''                # an empty annotation
        else annotations[i-1] <- intToUtf8 (tal[(endings[i-1]+1):(endings[i]-1)])
    }
    list (onset=onset, duration=duration, annotations=annotations)
}

doMergeASignals <- function (ASignals, annoSN) {
    for (sn in 1:length(ASignals)) {
        ASignals[[sn]]$annotations$fromSignal <- annoSN[sn]
    }
    ASignal1 <- ASignals[[1]]
    mergedAnnos <- ASignal1$annotations
    l <- length(ASignals)
    if (l>1) for (i in 2:l) {
        mergedAnnos <- rbind (mergedAnnos, ASignals[[i]]$annotations)
    }
    ASignal1$annotations    <- mergedAnnos
    ASignal1$signalNumber   <- annoSN
    ASignal1
}

# ------------------------------------------------------------------------------
#                          D signal functions
# ------------------------------------------------------------------------------
fragmentsToSignal <- function (fsignal, sSelFromRRT, sSelFromRS, sSelTillRS) {
    fsignal$signal      <- fragmentsSignalToSignal (fsignal$rFragments, sSelFromRS, sSelTillRS)
    if (length(fsignal$signal)) {
        fsignal$start       <- sSelFromRRT
        fsignal$fromSample  <- sSelFromRS
    }
    else {
        fsignal$start       <- as.numeric (NA)
        fsignal$fromSample  <- as.numeric (NA)
    }
    class(fsignal)  <- 'ebdfCSignal'
    fsignal
}

fragmentsSignalToSignal <- function (rFragments, sSelFromRS, sSelTillRS) {
    signal      <- integer (length=0)
    nSamples    <- max (0, sSelTillRS - sSelFromRS + 1)
    if (nSamples) {
        signal          <- rep (as.numeric(NA), nSamples)                       # memory check ??
        if (length(rFragments)) for (fn in 1:length(rFragments)) {
            rFragment   <- rFragments[[fn]]
            fStartRS    <- rFragment$fsFromRS
            fFromFS     <- fStartRS - sSelFromRS + 1                                # from in from-till signal fragment
            fTillFS     <- fFromFS + length (rFragment$signal) - 1
            signal[fFromFS:fTillFS] <- rFragment$signal
        }
    }
    signal
}

concatenateFragments <- function (rFragments, recordL) {
    nOld <- length (rFragments)
    fragments <- vector("list", length=0)
    if (nOld) {
        oldToNew    <- getOldToNew (rFragments, recordL)
        fragments   <- getNewFragments (rFragments, oldToNew)
    }
    fragments
}

getOldToNew <- function (rFragments, recordL) {
    nOld        <- length (rFragments)
    oldToNew    <- integer (length=nOld)
    tol         <- 1e-6

    newCnt      <- 1
    oldToNew[1] <- newCnt
    if (nOld > 1) for (oldCnt in 2:nOld) {
        s1 <- rFragments[[oldCnt-1]]$drsFromRS + recordL
        s2 <- rFragments[[oldCnt  ]]$drsFromRS
        if (s1 != s2) newCnt <- newCnt +1
        oldToNew [oldCnt] <- newCnt
    }
    oldToNew
}

getNewFragments <- function (rFragments, oldToNew) {                            # per signal
    nOld        <- length (oldToNew)
    nNew        <- oldToNew [nOld]

    oldFLength  <- sapply (rFragments, function (x){length(x$signal)} )
    newFragments<- vector("list", length=nNew)                                  # initialise for new fragments
    # class(newFragments) <- 'ebdfFragments'

    for (nfi in 1:nNew) {
        toCopyL <- oldToNew == nfi
        toCopyN <- which (toCopyL)

        newSignal   <- numeric (length = sum(oldFLength[toCopyL]))
        newFrom     <- 1
        for (j in toCopyN) {
            nextFrom <- newFrom + oldFLength[j]
            newSignal [newFrom:(nextFrom-1)] <- rFragments[[j]]$signal
            newFrom <- nextFrom
        }

        fromART <- rFragments[[toCopyN[1]]]$fsFromART
        fromRS  <- rFragments[[toCopyN[1]]]$fsFromRS
        fromRRT <- rFragments[[toCopyN[1]]]$fsFromRRT
        newFragments[[nfi]] <- list (start=fromART, fromSample=fromRS,
                                     recordingStart=fromRRT, signal= newSignal)
    }
    newFragments
}

