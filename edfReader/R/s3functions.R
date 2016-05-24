#                           S3 function for edfReader
#
# Purpose   :   print and summarize functions for ebdf objects
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
#
# History    :
#   Oct15 - Created
#   Mar16 - Revised, version 1.1.0
#   Apr16 - version 1.1.1, no changes
# ------------------------------------------------------------------------------
#                           s3 header functions
#                        objects: ebdfHeader & ebdfSHeader
# ------------------------------------------------------------------------------
#' @export
print.ebdfHeader <- function (x, ...) {                                         # excluding edf/bdf encoding data
    cat (' Patient               :', x$patient, '\n')
    cat (' RecordingId           :', x$recordingId, '\n')
    cat (' StartTime             :', posixltToChar(x$startTime), '\n')
    cat (' Continuous recording  :', x$isContinuous, '\n')
    labels <- x$sHeaders$label
    cat (' Signal labels         :', paste(labels, collapse = ', '), '\n')
    names <- row.names(x$sHeaders)
    if (sum(labels != names)) {                                                 # different names, i.e. non-unique labels
        cat (' R signal names        :', paste(names, collapse = ', '), '\n')
        cat (' Note                  : duplicate label names')
    }
}

#' @export
summary.ebdfHeader <- function (object, ...) {
    suffix <- ifelse (object$isPlus, '+', '')
    ft     <- paste (object$fileType, suffix, sep='')
    cat (' File name             :', object$fileName, '\n')
    cat (' File type             :', ft, '\n')
    cat (' Version               :', object$version, '\n')
    cat (' Patient               :', object$patient, '\n')
    cat (' RecordingId           :', object$recordingId, '\n')
    cat (' StartTime             :', posixltToChar (object$startTime), '\n')
    cat (' Continuous recording  :', object$isContinuous, '\n')
    if (!is.na(object$recordedPeriod)) {
        cat (' Recorded period       :', secPlusHHMMSS (object$recordedPeriod), '\n')
    }
    aCh <- sum (object$sHeaders$isAnnotation)
    cat (' Ordinary signals      :', object$nSignals - aCh, '\n')
    cat (' Annotation signals    :', aCh, '\n')
    labels <- object$sHeaders$label
    cat (' Signal labels         :', paste(labels, collapse = ', '), '\n')
    names <- row.names(object$sHeaders)
    if (sum(labels != names)) {                                                 # different names, i.e. non-unique labels
        cat (' R signal names        :', paste(names, collapse = ', '), '\n')
        cat (' Note                  : duplicate label names')
    }
}

#' @export
print.ebdfSHeaders <- function (x, ...) {
    labels <- x$label
    cat (' Signal labels         :', paste(labels, collapse = ', '), '\n')
    names <- row.names(x)
    if (sum(labels != names)) {                                                 # different names, i.e. non-unique labels
        cat (' R signal names        :', paste(names, collapse = ', '), '\n')
        cat (' Note                  : duplicate label names')
    }
}

#' @export
summary.ebdfSHeaders <- function (object, ...) {
    sns     <- seq (to=nrow(object))
    labels  <- object$label
    names   <- row.names(object)
    if (!sum(labels == names)) {
        sdf <- data.frame (signal=sns, label=labels, transducer=object$transducerType,
                           sampleRate=object$sRate, preFilter= object$preFilter, stringsAsFactors=FALSE)
    } else {
        sdf <- data.frame (signal=sns, label=labels, name=names, transducer=object$transducerType,
                           sampleRate=object$sRate, preFilter= object$preFilter, stringsAsFactors=FALSE)
    }
    annons <- which(sdf$label=='EDF Annotations' | sdf$label=='BDF Annotations')
    if (length(annons)) {
        sdf$transducer[annons]  <- ' '
        sdf$sampleRate[annons]  <- NA
        sdf$preFilter[annons]   <- ' '
    }
    sdf
}
# ------------------------------------------------------------------------------
#                            s3 lists of signals
#                            object: ebdfSignals
# ------------------------------------------------------------------------------
#' @export
print.ebdfSignals <- function (x, ...) {
    printSummarySignals (x, isSumm=FALSE)
}

#' @export
summary.ebdfSignals <- function (object, ...) {
    printSummarySignals (object, isSumm=TRUE)
}

printSummarySignals <- function (signals, isSumm) {
    asn  <- logical( length = length(signals))
    csn  <- asn
    fsn  <- asn                     # +D read with frogmented=TRUE
    for (sn in 1:length(signals)) {
        if (class (signals[[sn]]) == 'ebdfASignal') asn[sn] <- TRUE
        if (class (signals[[sn]]) == 'ebdfCSignal') csn[sn] <- TRUE
        if (class (signals[[sn]]) == 'ebdfFSignal') fsn[sn] <- TRUE
    }

    printSummaryCommonSignalsData      (signals, isSumm=isSumm)
    if (sum(asn)) printSummaryASignals (signals, asn=asn, isSumm=isSumm)
    if (sum(csn)) printSummaryCSignals (signals, csn=csn, isSumm=isSumm)
    if (sum(fsn)) printSummaryFSignals (signals, fsn=fsn, isSumm=isSumm)
}

printSummaryCommonSignalsData <- function (signals, isSumm) {
    # osn1, asn1, fns1 : the first in signals if present else 0
    if (isSumm) {
        cat (' StartTime               :', posixltToChar(signals[[1]]$startTime), '\n')
    }
}

printSummaryASignals <- function (signals, asn,  isSumm) {
    nr  <- sum(asn)
    cat (ifelse(nr==1,'Annotation signal:', 'Annotation signals:' ), '\n')
    sns     <- character (length=nr)
    from    <- character (length=nr)
    till    <- character (length=nr)
    n       <- integer   (length=nr)
    r <- 1
    for (sn in 1:length(signals)) if (asn[sn]) {
        sns[r]      <- paste (signals[[sn]]$signalNumber, collapse=',')
        fromTillN   <- getAnnotsSumASignal (signals[[sn]])
        from[r]     <- secPlusHHMMSS (fromTillN [1])
        till[r]     <- secPlusHHMMSS (fromTillN [2])
        n[r]        <- fromTillN [3]
        r           <- r +1
    }
    pdf <- data.frame (signal=sns, annotations=n , from=from, till=till, stringsAsFactors = FALSE)
    print.data.frame  (pdf)
}

printSummaryCSignals <- function (signals, csn, isSumm) {
    cat (ifelse(sum(csn)==1,'Ordinary signal:', 'Ordinary signals:' ), '\n')
    printSummaryOSignalCommon (signals[[which.max(csn)]], isSumm, commonsOnly=TRUE)
    printSummaryOSignalsList (signals, csn)
}

printSummaryFSignals <- function (signals, fsn, isSumm) {
    cat (ifelse(sum(fsn)==1,'Ordinary signal:', 'Ordinary signals:' ), '\n')
    oSignal <-signals[[which.max(fsn)]]
    printSummaryOSignalCommon (oSignal, isSumm, commonsOnly=TRUE)
    printSummaryOSignalsList (signals, fsn)
}

printSummaryOSignalsList <- function (signals, osn) {
    nr      <- sum(osn)
    oLabels <- sapply (signals, function (x){x$label})[osn]
    oNames  <- names(signals)[osn]
    dubs    <- sum (oLabels != oNames)

    sns  <- integer   (length=nr)
    tds  <- character (length=nr)
    srs  <- numeric   (length=nr)
    sms  <- integer   (length=nr)
    pfs  <- character (length=nr)
    fRS  <- integer   (length=nr)
    fRT  <- numeric   (length=nr)
    r <- 1
    for (sn in 1:length(signals)) if (osn[sn]) {
        sns[r] <- signals[[sn]]$signalNumber
        tds[r] <- signals[[sn]]$transducer
        srs[r] <- signals[[sn]]$sRate
        sms[r] <- length(signals[[sn]]$signal)
        pfs[r] <- signals[[sn]]$preFilter
        r <- r +1
    }
    if (!dubs) {
        psdf <- data.frame (signal=sns, label=oLabels,
                            transducer=tds, sampleRate=srs, samples=sms,
                            preFilter=pfs,stringsAsFactors = FALSE)
    } else {
        psdf <- data.frame (signal=sns, label=oLabels, name=oNames,
                            transducer=tds, sampleRate=srs, samples=sms,
                            preFilter=pfs,stringsAsFactors = FALSE)
    }
    print.data.frame  (psdf)
}

# ------------------------------------------------------------------------------
#                           s3 signal functions
#            objects: ebdfASignal, ebdfCSignal, and ebdfFSignal
# ------------------------------------------------------------------------------
#' @export
print.ebdfASignal <- function (x, ...) {
    printSummaryASignal (aSignal=x, isSumm=TRUE)
}

#' @export
summary.ebdfASignal <- function (object, ...) {
    printSummaryASignal (aSignal=object, isSumm=TRUE)
}

printSummaryASignal <- function (aSignal, isSumm) {
    cat (" Signal number           :", aSignal$signalNumber, '\n')
    cat (" Label                   :", aSignal$label, '\n')
    cat (' StartTime               :', posixltToChar(aSignal$startTime), '\n')


    rss <- sum (aSignal$annotations$isRecordStart)
    an  <- nrow (aSignal$annotations) - rss
    cat (" Record start specs      :", rss, '\n')
    cat (" Annotations             :", an,  '\n')

    fromTillN   <- getAnnotsSumASignal (aSignal)
    cat (" Time first annotation   :", secPlusHHMMSS (fromTillN [1]), '\n')
    cat (" Time last annotation    :", secPlusHHMMSS (fromTillN [2]), '\n')
}

#' @export
print.ebdfCSignal <- function (x, ...) {
    printSummaryOSignalCommon (oSignal=x, isSumm=TRUE, commonsOnly=FALSE)
}

#' @export
summary.ebdfCSignal <- function (object, ...) {
    printSummaryOSignalCommon (oSignal=object, isSumm=TRUE, commonsOnly=FALSE)
    cat (" Signal summary          :\n")
    summary (object$signal)
}

#' @export
print.ebdfFSignal <- function (x, ...) {
    printSummaryOSignalCommon (oSignal=x, isSumm=TRUE, commonsOnly=FALSE)
}

#' @export
summary.ebdfFSignal <- function (object, ...) {
    printSummaryOSignalCommon (oSignal=object, isSumm=TRUE, commonsOnly=FALSE)
    nFragments <- length (object$fragments)
    rows <- ifelse(nFragments>10, 5, nFragments)                                # 5 if > 10
    psdf <- getFragmentSummaries (object, rows)
    print (psdf)
    if (nFragments > rows) cat ("....", nFragments-rows, "more\n")
    cat ("All fragments:\n")
    print (getFragmentsSummary(object))
}

printSummaryOSignalCommon <- function (oSignal, isSumm, commonsOnly) {

    wRecTillRRT <- oSignal$totalPeriod
    readFromRRT <- max (0,           oSignal$from)
    readTillRRT <- min (wRecTillRRT, oSignal$till)

    if (!commonsOnly) {
        cat (" Signal number           :", oSignal$signalNumber, '\n')
        cat (" Label                   :", oSignal$label, '\n')
    }

    if (!commonsOnly) {                                                         # exception
        cat (' StartTime               :', posixltToChar(oSignal$startTime), '\n')
    }

    cat (" Continuous recording    :", oSignal$isContinuous,   '\n')
    if (!isSumm) {
        pReadTxt    <- getPeriodFromText (fromRRT = readFromRRT, tillRRT = readTillRRT,
                                          wRectillRRT = wRecTillRRT )
        cat (" Period read             :", pReadTxt, '\n')
    } else {                                                                    # a summary
        if (!oSignal$isContinuous) {
            cat (" Total period            :", secPlusHHMMSS (wRecTillRRT), '\n' )
            pReadLabel <- " Read from total period  :"
        } else {
            pReadLabel <- " Period read             :"
        }
        cat (' Recorded period         :', secPlusHHMMSS (oSignal$recordedPeriod), '\n')
        pReadTxt    <- getPeriodFromText (fromRRT = readFromRRT, tillRRT = readTillRRT)
        cat (pReadLabel, pReadTxt, '\n')

        if (!commonsOnly) {
            cat (" Transducer              :", oSignal$transducerType, '\n')
            cat (" Range                   :", oSignal$range, '\n')
            cat (" Prefilter               :", oSignal$preFilter, '\n')
        }

        cat (" Bits per sample         :", oSignal$sampleBits, '\n')
        if (!commonsOnly) {
            cat (" Sample rate             :", oSignal$sRate, '\n')
        }

        if (class (oSignal) == 'ebdfFSignal') {
            if (commonsOnly) {
                fragmentsLabel <- " Fragments per signal    :"
            } else {
                fragmentsLabel <- " Number of fragments     :"
            }
            cat (fragmentsLabel, length (oSignal$fragments), '\n')
        } else {
            cat (" Number of samples       :", length (oSignal$signal), '\n')
        }
    }
}

getFragmentsSummary <- function (object) {
    nFragments  <- length (object$fragments)
    totalLength <- 0
    for (i in 1:nFragments) {
        totalLength <- totalLength + length(object$fragments[[i]]$signal)
    }
    totalSignal <- numeric (length = totalLength)
    from <- 1
    for (i in 1:nFragments) {
        l <- length(object$fragments[[i]]$signal)
        totalSignal[from:(from+l-1)] <- object$fragments[[i]]$signal
        from <- from + l
    }
    summary (totalSignal)
}

getFragmentSummaries <- function (object, rows) {
    n       <- integer(length=rows)
    starts  <- numeric(length=rows)
    samples <- integer(length=rows)
    mins    <- numeric(length=rows)
    q1s     <- numeric(length=rows)
    medians <- numeric(length=rows)
    means   <- numeric(length=rows)
    q3s     <- numeric(length=rows)
    maxs    <- numeric(length=rows)
    for (i in 1:rows) {
        n[i]        <- i
        starts[i]   <- object$fragments[[i]]$start
        samples[i]  <- length(object$fragments[[i]]$signal)
        summ        <- summary (object$fragments[[i]]$signal)
        mins[i]     <- summ ["Min."]
        q1s[i]      <- summ ["1st Qu."]
        medians[i]  <- summ ["Median"]
        means[i]    <- summ ["Mean"]
        q3s[i]      <- summ ["3rd Qu."]
        maxs[i]     <- summ ["Max."]
    }
    psdf <- data.frame (fragment=n, start=starts, samples=samples, "Min."=mins, "1st Qu."=q1s,
                        "Median"=medians, "Mean"=means, "3rd Qu."=q3s, "Max."=maxs, stringsAsFactors = FALSE)
    psdf
}

# ------------------------------------------------------------------------------
#                        common functions
# ------------------------------------------------------------------------------

# period read :  xxxx sec from yyy = recorded period / from start / till end of recording
getPeriodFromText <- function (fromRRT, fromRS=NULL, tillRRT, wRectillRRT=NULL) {

    periodTxt <- secPlusHHMMSS (tillRRT - fromRRT)
    if (fromRRT) {
        periodTxt <- paste (periodTxt, ' starting at ', secPlusHHMMSS (fromRRT), sep='')
    }
    if (!is.null(fromRS)) {
        periodTxt <- paste (periodTxt, ' at  sample ', fromRS, sep='')
    }

    if ( !is.null (wRectillRRT)) {                                              # compose additional text
        isFromStart <- fromRRT == 0
        isTillEnd   <- wRectillRRT <= tillRRT
        if (isFromStart & isTillEnd) {                                          # whole recoding
            periodTxt <- paste (periodTxt, ' = whole recording', sep='')
        } else if (isFromStart & !isTillEnd ) {                                 # from start till
            periodTxt <- paste (periodTxt, ' =  from start', sep='')
        } else if (!isFromStart & isTillEnd ) {                                 # from tlll end
            periodTxt <- paste (periodTxt, ' till the end', sep='')
        }
    }
    periodTxt
}

getAnnotsSumASignal <- function (aSignal) {
    annots  <- aSignal$annotations
    onsets  <- annots [!annots$isRecordStart, 'onset']
    n       <- length (onsets)
    from    <- onsets [1]
    till    <- onsets [n]
    c (from, till, n)
}

secPlusHHMMSS <- function (sec) {
    note <- ''
    if (sec > 60) {
        hhmmss <- secToHHMMSS (sec)
        note  <- paste (' (= ', hhmmss, ')', sep='')
    }
    paste (sec, ' sec' ,note, sep='')
}

secToHHMMSS <- function (sec) {
    hhmmss <- ''
    if (sec) {
        m   <- sec %/% 60    # rounded to minutes
        hhmmss <- sprintf ('%02d:%02d:%09.6f', m%/%60, m%%60, sec%%60)
        # remove trailing '0's and '.'
        n <- nchar(hhmmss)
        if (n > 8) {
            while (substr(hhmmss, n, n) == '0') n <- n-1
            if    (substr(hhmmss, n, n) == '.') n <- n-1
            hhmmss <- substr(hhmmss, 1, n)
        }
    }
    hhmmss
}

posixltToChar <- function (t) {
    dt <- format (t, format="%Y-%m-%d %H:%M:%OS9",  usetz = FALSE)
    n <- nchar(dt)
    # remove trailing '0's and '.'
    if (n > 19) {
        while (substr(dt, n, n) == '0') n <- n-1
        if    (substr(dt, n, n) == '.') n <- n-1
        dt <- substr(dt, 1, n)
    }
    dt
}
