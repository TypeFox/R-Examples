#
# Purpose   :   Test the reading of complete EDF+D and BDF+D files with edfReader
#
#               The +D file should be derived from an +C file by ommitting data records
#               The +D file is tested againts the original +C file
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
#   Mar16 - Created for version 1.1.0
# ------------------------------------------------------------------------------

require (testthat)
require (edfReader)
context ("Compare reading whole files with the export from EDFBrowser.")

libDir <- paste (system.file("extdata", package="edfReader"), '/', sep='')
DFns <- 'edfPlusD.edf'                           # a subset of edfPlusC
# the source, i.e. the original +C files in the same order:
SFns <- 'edfPlusC.edf'

# edfPlusC.edf is a truncated version of a test_generator_2 test file
# the original file can be found at: http://www.teuniz.net/edf_bdf_testfiles
# edfPlusD.edf is derived from edfPlusC.edf by deleting some data records

fnsN <- length(DFns)

DFFns <- character (length = fnsN)
for (i in 1:fnsN) DFFns[i] <- paste (libDir, DFns[i], sep='')
SFFns <- character (length = fnsN)
for (i in 1:fnsN) SFFns[i] <- paste (libDir, SFns[i], sep='')

DHdrs <- vector (mode='list', length = fnsN)
for (i in 1:fnsN) DHdrs[[i]] <- readEdfHeader(DFFns[i])
SHdrs <- vector (mode='list', length = fnsN)
for (i in 1:fnsN) SHdrs[[i]] <- readEdfHeader(SFFns[i])

#                          test header and sheader
# ------------------------------------------------------------------------------
testDHeader <- function (SHdr, DHdr) {
    # apply changes to SHdr and test equallity with DHdr
    NHdr <- SHdr
    NHdr$fileName       <- DHdr$fileName

    if (NHdr$fileType == 'EDF') NHdr$reserved   <- "EDF+D"
    else                        NHdr$reserved   <- "BDF+D"
    NHdr$isContinuous       <- FALSE
    NHdr$nRecords           <- DHdr$nRecords
    NHdr$recordedPeriod     <- DHdr$recordedPeriod
    NHdr$sHeaders$sLength   <- NHdr$sHeaders$samplesPerRecord * NHdr$nRecords
    test_that ("DHdr is equal to NHdr", {
        expect_equal(NHdr, DHdr)
    })
}
#                          test signals
# ------------------------------------------------------------------------------

testDOSignals  <- function (DHdr, SSignals, DSignalsC, DSignalsF, SFn) {
    nSignals    <- length (SSignals)
    SOSignalsL  <- sapply (SSignals,  function(X){!X$isAnnotation})
    DOSignalsCL <- sapply (DSignalsC, function(X){!X$isAnnotation})
    DOSignalsFL <- sapply (DSignalsF, function(X){!X$isAnnotation})

    that <- paste ("DOS1: ", SFn, ": Same signals: number, names, ordinary", sep='')
    test_that (that, {
        expect_equal (length(SSignals), length(DSignalsC))
        expect_equal (length(SSignals), length(DSignalsF))
        expect_equal (names (SSignals), names (DSignalsC))
        expect_equal (names (SSignals), names (DSignalsF))
        expect_equal (SOSignalsL, DOSignalsCL)
        expect_equal (SOSignalsL, DOSignalsFL)
    })

    that <- paste ("DOS2: ", SFn, ": Total D signals length", sep='')
    test_that (that, {
        for (sn in 1: nSignals) if (SOSignalsL[sn]) {
            rns     <- names (SSignals)[sn]
            lS      <- DHdr$sHeaders[rns, 'sLength']
            notNa   <- !is.na (DSignalsC[[sn]]$signal)
            lC      <- sum (notNa)
            expect_equal(lS, lC)

            fragments <- DSignalsF[[sn]]$fragments
            lF <- 0
            for (fi in 1:length(fragments)) {
                lF <- lF + length(fragments[[fi]]$signal)
            }
            expect_equal(lS, lF)
        }
    })

    that <- paste ("DOS3: ", SFn, ": D signals values", sep='')
    test_that (that, {
        for (sn in 1: nSignals) if (SOSignalsL[sn]) {
            notNa   <- !is.na (DSignalsC[[sn]]$signal)
            expect_equal (DSignalsC[[sn]]$signal[notNa], SSignals[[sn]]$signal[notNa])
            SSignal <- SSignals[[sn]]$signal
            for (i in 1:length(SSignal)) {
                SSignal[i] <- ifelse(notNa[i], SSignal[i], NA)
            }
            expect_equal (DSignalsC[[sn]]$signal, SSignal)

            fragments <- DSignalsF[[sn]]$fragments
            fragmentsG <<- fragments
            for (fi in 1:length(fragments)) {
                fromRS <- fragments[[fi]]$fromSample
                tillRS <- fromRS + length (fragments[[fi]]$signal) - 1
                info   <- paste ('fragment=', fi, sep='')
                expect_equal (fragments[[fi]]$signal, SSignals[[sn]]$signal[fromRS:tillRS], info=info)
            }
        }
    })
}

testDAnnotations <- function (SHdr, DHdr, SSignals, DSignalsC, DSignalsF, SFn) {
    nSignals    <- length (SSignals)
    SASignalsL  <- sapply (SSignals,  function(X){X$isAnnotation})
    DASignalsCL <- sapply (DSignalsC, function(X){X$isAnnotation})
    DASignalsFL <- sapply (DSignalsF, function(X){X$isAnnotation})

    that <- paste ("DAS1: ", SFn, ": Same signals: number, names, ordinary", sep='')
    test_that (that, {
        expect_equal (length(SSignals), length(DSignalsC))
        expect_equal (length(SSignals), length(DSignalsF))
        expect_equal (names (SSignals), names (DSignalsC))
        expect_equal (names (SSignals), names (DSignalsF))
        expect_equal (SASignalsL, DASignalsCL)
        expect_equal (SASignalsL, DASignalsFL)
        expect_equal (sum (SASignalsL), 1)                                      # merging assumed
    })

    that <- paste ("DAS2: ", SFn, ": Fragmented/nonFragmented: same annotations", sep='')
    test_that (that, {
        for (sn in 1: nSignals) if (SASignalsL[sn]) {
            expect_equal (DSignalsC[[sn]]$annotations,DSignalsF[[sn]]$annotations)
        }
    })

    getRecordsIncludedinD <- function () {
        os1n        <- which.min(SASignalsL)
        os1Name     <- names(SSignals)[os1n]
        sPerRec     <- DHdr$sHeaders[os1Name,]$samplesPerRecord
        includedL   <- logical (length = SHdr$nRecords)
        fragments   <- DSignalsF[[os1n]]$fragments
        for (fn in 1:length(fragments)) {
            fragment <- fragments[[fn]]
            fromRec  <- (fragment$fromSample -1) / sPerRec + 1
            fRecsN   <- length (fragment$signal) / sPerRec
            for (i in 1:fRecsN) {
                includedL[fromRec+i-1] <- TRUE
            }
        }
        includedL
    }

    testSignalAnnotations <- function (sn) {
        SASignal    <- SSignals [[sn]]
        DASignal    <- DSignalsF[[sn]]
        includedL   <- getRecordsIncludedinD()
        fromDataRec <- which.max(includedL)
        tillDataRec <- length(includedL) - which.max(rev(includedL)) + 1

        DFromT      <- DHdr$recordDuration * (fromDataRec -1)
        test_that("DAS3: From is equal or correctly adusted", {
            expect_equal (SASignal$from + DFromT, DASignal$from)
        })
        test_that("DAS4 Till is equal", {
            expect_equal (SASignal$till , DASignal$till)
        })

        SARecs      <- SASignal$annotations$record
        SARowsInDL  <- logical(length = length(SARecs))
        for (rn in 1:length(SARecs)) SARowsInDL[rn] <- includedL[SARecs[rn]]

        if (sum (SARowsInDL) == 0)  {
            test_that("DAS5: No annotations included", {
                expect_equal (nrow(DASignal$annotations), 0)
            })
        } else {
            SAnnots <- SASignal$annotations[SARowsInDL,]
            SAnnots$record <- NULL
            DAnnots <- DASignal$annotations
            DAnnots$record <- NULL
            test_that("DAS6: Annotations included are equal", {
                expect_equivalent (SAnnots, DAnnots)
            })
        }
    }

    for (sn in 1: nSignals) if (SASignalsL[sn]) {
        testSignalAnnotations (sn)
    }
}
#                              test D files
# ------------------------------------------------------------------------------
testDFile <- function (n) {
    cat (DFns[n], ": Testing header\n")
    SHdr <- SHdrs[[n]]
    DHdr <- DHdrs[[n]]
    testDHeader   (SHdr = SHdr, DHdr = DHdr)

    cat (DFns[n], ": Testing signals\n")
    SSignals    <- readEdfSignals (SHdr, simplify=FALSE)
    DSignalsC   <- readEdfSignals (DHdr, fragments=FALSE, simplify=FALSE)
    DSignalsF   <- readEdfSignals (DHdr, fragments=TRUE,  simplify=FALSE)

    testDOSignals (DHdr=DHdr, SSignals=SSignals, DSignalsC=DSignalsC, DSignalsF=DSignalsF, SFn=SFns[n])
    testDAnnotations (SHdr=SHdr, DHdr=DHdr, SSignals=SSignals, DSignalsC=DSignalsC, DSignalsF=DSignalsF, SFn=SFns[n])
}

testAllDFiles <- function () {
    testDFile (1)
}

testAllDFiles()
