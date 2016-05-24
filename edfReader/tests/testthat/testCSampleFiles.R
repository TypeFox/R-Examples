#
# Purpose   :   Test the reading of complete continuous .edf(+)/.bdf(+) files with edfReader
#
#               The readings are compared with the imported ASCII export
#               produced with EDFBrowser (by Teunis van Beelen)
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
#               along with edf package for R.  If not, see <http://www.gnu.org/licenses/>.
#
# History    :
#   Jan16 - Created
#   Mar16 - revised for verion 1.1.0;  test for +D file moved to another file;
#           suppport for first record onset != 0
# ------------------------------------------------------------------------------

require (testthat)
require (edfReader)
context ("Compare reading whole files with the export from EDFBrowser.")

libDir <- paste (system.file("extdata", package="edfReader"), '/', sep='')
sFns <- c('edfPlusC.edf',
          'bdfPlusC.bdf',          # same signals as in edfPlusC.edf, but 24 bits
          'edfAnnonC.edf')         # with 2 annotation signals, to test annotations only

# edfPlusC.edf and bdfPlusC.bdf are truncated versions of the test_generator_2 test file
# the original files can be found at: http://www.teuniz.net/edf_bdf_testfiles
# edfAnnonC.edf  is a truncated version of the test_generator8.edf test file
# test_generator8.edf has been received from Teunis van Beelen via private communications

sFFns <- character (length = length(sFns))
for (i in 1:length(sFns)) sFFns[i] <- paste (libDir, sFns[i], sep='')

sHdrs <- vector (mode='list', length = length(sFns))
for (i in 1:length(sFns)) sHdrs[[i]] <- readEdfHeader(sFFns[i])

#                          importBrowserExport
# ------------------------------------------------------------------------------
importBrowserExport <- function (sFn) {
    fn      <- substring (sFn, 1, nchar(sFn)-4)
    hdrFFn  <- paste (libDir,fn, '_header.txt',      sep='')
    sHdrFFn <- paste (libDir,fn, '_signals.txt',     sep='')
    annoFFn <- paste (libDir,fn, '_annotations.txt', sep='')
    dataFFn <- paste (libDir,fn, '_data.txt',        sep='')

    hdr     <- read.csv(hdrFFn,  stringsAsFactors = FALSE)
    hdr[is.na(hdr)] <-""                                                        # na.strings = "" in read.csv didn't seem to do thejob
    sHdr     <- read.csv(sHdrFFn, stringsAsFactors = FALSE)
    sHdr[is.na(sHdr)] <-""                                                      # na.strings = "" in read.csv didn't seem to do thejob
    if (file.exists(annoFFn)) {
        anno    <- read.csv(annoFFn, stringsAsFactors = FALSE)
        anno    <- anno[order(anno$Onset),]
        anno$Duration <- as.numeric(anno$Duration)                              # none specified, duration will be imported as logical
    } else {
        anno <- NULL
    }

    if (file.exists(dataFFn)) {
        data     <- read.csv(dataFFn, stringsAsFactors = FALSE)
    } else {
        data <- NULL
    }
    list (header = hdr, sHeaders=sHdr, annotations=anno, signals=data)
}
#                          test header and sheader
# ------------------------------------------------------------------------------
testCHeader <- function (edfHdr, expHdr, sFn) {
    edfVersion <- edfHdr$version                                                # version
    oSignalsL <- !edfHdr$sHeaders$isAnnotation
    nOSignals <- sum (!edfHdr$sHeaders$isAnnotation)
    if (edfVersion == "\"255\"BIOSEMI") edfVersion <- ".BIOSEMI"
    that <- paste ("H1: ", sFn, ": Equal file headers (ex signal info)", sep='')
    test_that(that, {
        expect_equal (edfVersion, as.character (expHdr$Version))
        expect_equal (edfHdr$patient, trimws(expHdr$Subject))                   # patient
        expect_equal (edfHdr$recordingId, trimws(expHdr$Recording))             # recording id
        expect_equal (format (edfHdr$startTime, "%d.%m.%y"), expHdr$Startdate)  # start date
        expect_equal (format (edfHdr$startTime, "%H.%M.%S"), expHdr$Startime)   # start time
        expect_equal (edfHdr$headerLength, expHdr$Bytes)                        # header length in octets
        expect_equal (edfHdr$nRecords, expHdr$NumRec)                           # number of data records
        expect_equal (edfHdr$recordDuration, expHdr$Duration)                   # record duration
        if (edfHdr$isPlus) {
            expect_equal (nOSignals, expHdr$NumSig)                             # number of signals
        } else {
            expect_equal (edfHdr$nSignals, expHdr$NumSig)                       # number of signals
        }
    })
}

testCSHeaders <- function (edfHdr, edfSHdr, expSHdr, sFn) {
    # testSHeaders (hdr, hdr$sHeaders, exp$sHeaders, sFn=sFns[n])
    labels <- edfSHdr$label
    names <- row.names (edfSHdr)
    diffs <- sum (labels!=names)
    that <- paste ("SH1: ", sFn, ": Name-label equallity", sep='')
    if (sFn==sFns[1] | sFn==sFns[2]) {                                          # sFns[3] contains non-unique labels
        expect_equal (diffs, 0)
    }
    # EDFBrowser doesn't export sHeader data for annotation signals in EDF files
    # So, delete these rows form edfSHdr, if exist
    if ( edfHdr$isPlus & sum (edfSHdr$isAnnotation)) {
        edfSHdr <- edfSHdr [!edfSHdr$isAnnotation, ]
    }

    # "," is the export field separator, therfore EDFBrowser changed it to"'"
    # So, apply the same change as EDFBrowser
    edfTranducer <- gsub(",", "'", edfSHdr$transducerType)

    that <- paste ("SH2: ", sFn, ": Equal signal headers", sep='')
    test_that (that, {
        expect_equal (nrow(edfSHdr), nrow(expSHdr))
        if (nrow(edfSHdr)) {
            expect_equal (edfSHdr$label, trimws(expSHdr$Label))                 # label
            expect_equal (edfTranducer, trimws(expSHdr$Transducer))             # Transducer type
            expect_equal (edfSHdr$physicalDim, trimws(expSHdr$Units))           # Dimension
            expect_equal (edfSHdr$physicalMin, expSHdr$Min)                     # Physical min
            expect_equal (edfSHdr$physicalMax, expSHdr$Max)                     # Physical max
            expect_equal (edfSHdr$digitalMin, expSHdr$Dmin)                     # Digital min
            expect_equal (edfSHdr$digitalMax, expSHdr$Dmax)                     # Digital max
            expect_equal (edfSHdr$preFilter, trimws(expSHdr$PreFilter))         # pre filters
            expect_equal (edfSHdr$samplesPerRecord, expSHdr$Smp.Rec)            # Samples per record
            expect_equal (edfSHdr$reserved, trimws(expSHdr$Reserved))           # Reserved
        }
    })
}
#                          test signals
# ------------------------------------------------------------------------------
testCSignals <- function (edfHdr, edfSignals, expSignals, sFn) {
    cSignalsL   <- sapply (edfSignals, function(X){!X$isAnnotation})
    nCSignals   <- sum(cSignalsL)
    cStoS       <- which (cSignalsL)   # ordinary Signal to Signal

    nExpSignals <- length(expSignals)-1

    if (nCSignals) {
        that <- paste ("COS1: ", sFn, ": Equal number of ordinary signals", sep='')
        test_that (that, {
            expect_equal (nCSignals, nExpSignals)
        })
        eSignals <- vector (mode='list', length = nCSignals)
        for (sn in 1:nCSignals) {
            # skip first colum in expSignals (contains the time)
            # remove 'na' from ASCII / data frame colomn due to
            #   -  different sample rates and/or
            #   -  intermediate time only rows (for non-continuous signals)
            s <- expSignals[,(sn+1)]
            s <- s[!is.na(s)]
            eSignals[[sn]] <- s                                                 # exp signals to compare with
        }
        for (csn in 1:nCSignals) {
            sn      <- cStoS[csn]
            cSigAll <- edfSignals[[sn]]$signal
            # cat ("C2: length(cSigAll)=", length(cSigAll), "sum(is.na(cSigAll)=", sum(is.na(cSigAll)) ,'\n')
            cSig    <- cSigAll[!is.na(cSigAll)]

            that <- paste ("COS2: ", sFn, ": Equal length for continuous signal ", sn, sep='')
            test_that (that, {
                expect_equal (length(cSig), length(eSignals[[csn]]))
            })

            that <- paste ("COS3: ", sFn, ": Signal length equal to length in signal header", sn, sep='')
            test_that (that, {
                expect_equal (length(cSig), edfHdr$sHeaders$sLength[1])
            })

            that <- paste ("COS4: ", sFn, ": recording period equal to recording period in header", sn, sep='')
            test_that (that, {
                expect_equal (edfSignals[[sn]]$recordedPeriod, edfHdr$recordedPeriod)
            })

            that <- paste ("COS5: ", sFn, ": (Almost) equal continuous signals for sn ", sn, sep='')
            test_that (that, {
                expect_equal (cSig, eSignals[[csn]], tolerance=5e-7)
            })
        }
    }
}

testCAnnotations <- function (edfHdr, edfSignals, expAnnos, mergeASignals=TRUE, sFn) {
    HAsignalsL  <- edfHdr$sHeaders$isAnnotation
    SAsignalsL  <- sapply (edfSignals, function(X) X$isAnnotation)
    HAsignalsN  <- sum(HAsignalsL)
    SAsignalsN  <- sum(SAsignalsL)
    if (HAsignalsN) {
        # Test signal numbers
        that <- paste ("CAS1: ", sFn, ": Annotation signal numbers", sep='')
        test_that(that, {
            expect_equivalent (which.max(SAsignalsL), which.max(HAsignalsL))
            if (mergeASignals) {
                expect_equal (SAsignalsN, 1)
            } else {
                expect_equal (SAsignalsL, HAsignalsL)
            }
        })
        aToS <- which (SAsignalsL)
        for (sn in 1:SAsignalsN) {
            annos <- edfSignals[[aToS[sn]]]
            that <- paste ("CAS2: ", sFn, " sn=", sn ,": Equal annotation signals", sep='')
            test_that (that, {
                expect_equal (nrow(expAnnos), nrow(annos$annotations))
                expect_equal (expAnnos$Onset, annos$annotations$onset)
                expect_equal (expAnnos$Duration, annos$annotations$duration)
                # EDFBrowser removes '"' characters
                expect_equal (expAnnos$Annotation, gsub('"', "", annos$annotations$annotation))
            })
        }
    }
}
#                              test files
# ------------------------------------------------------------------------------
testCFile <- function (n) {
    cat (sFns[n], ": Testing header\n")
    exp     <- importBrowserExport (sFns[n])
    hdr     <- readEdfHeader(sFFns[n])
    testCHeader   (hdr, exp$header, sFn=sFns[n])
    testCSHeaders (hdr, hdr$sHeaders, exp$sHeaders, sFn=sFns[n])

    cat (sFns[n], ": Testing signals\n")
    cSignals <- readEdfSignals (hdr, simplify=FALSE)
    if (!is.null(exp$signals)) {
        testCSignals    (hdr, cSignals, exp$signals, sFn=sFns[n])
    }
    if (!is.null(exp$annotations)) {
        testCAnnotations (hdr, cSignals, exp$annotations, mergeASignals=TRUE, sFn=sFns[n])
    }
}

testAllFiles <- function () {
    testCFile (1)
    testCFile (2)
    testCFile (3)
}

testAllFiles()
