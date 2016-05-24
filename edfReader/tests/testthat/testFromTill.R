#
# Purpose   :    Test the use of the 'from' anbd 'till' paraams in readEdfSignals()
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
#   Mar16 - revised for verion 1.1.0; support added for +D file & first record onset != 0
# ------------------------------------------------------------------------------
#
require (testthat)
require (edfReader)

context ("Compare reading whole files with the export from EDFBrowser.")

libDir <- paste (system.file("extdata", package="edfReader"), '/', sep='')

sFns <- c('edfPlusC.edf',
          'bdfPlusC.bdf',          # same signals as in edfPlusC.edf, but 24 bits
          'edfPlusD.edf',          # extracted from 'edfPlusC.edf'
          'edfAnnonC.edf')         # with 2 annotation signals, to test annotations only'

# edfPlusC.edf and bdfPlusC.bdf are truncated versions of the test_generator_2 test file
# the original files can be found at: http://www.teuniz.net/edf_bdf_testfiles
# edfAnnonC.edf  is a truncated version of the test_generator8.edf test file
# test_generator8.edf has been received from Teunis van Beelen via private communications

sFFns <- character (length = length(sFns))
for (i in 1:length(sFns)) sFFns[i] <- paste (libDir, sFns[i], sep='')

sHdrs <- vector (mode='list', length = length(sFns))
for (i in 1:length(sFns)) sHdrs[[i]] <- readEdfHeader(sFFns[i])

#                        test from-till for annotations
# ------------------------------------------------------------------------------
testASignalsFile <- function (fileNo) {
    hdr <- sHdrs[[fileNo]]
    cat ("testASignalsFile", sFns[fileNo], '\n')

    isAnnotation <- hdr$sHeaders$isAnnotation
    nASignals   <- sum(isAnnotation)
    if (!nASignals) return

    # select an annotation signal
    asn         <- sample.int(nASignals, 1)
    asn         <- which(isAnnotation)[asn]
    wholeASign  <- readEdfSignals(hdr, signals=asn)
    totalPeriod <- wholeASign$totalPeriod
    allAnnots   <- wholeASign$annotations
    allOnsets   <- allAnnots$onset
    nOnsets     <- length (allOnsets)

    test_that ("A1: [from,from) one annotation included", {
        from    <- allOnsets[1]
        till    <- allOnsets[1]          # will cause an error if equal to total recording time
        aSignal <- readEdfSignals (hdr, signals=asn, from=from, till=till, fragment=FALSE)
        aFSignal<- readEdfSignals (hdr, signals=asn, from=from, till=till, fragment=TRUE )
        expect_equal (aSignal, aFSignal)
        expect_equal (aSignal$from , from)
        expect_equal (aSignal$till , till)
        annotsInRange <- (from <= allOnsets) & (allOnsets <= till)
        expect_equal (aSignal$annotations, allAnnots[annotsInRange, ])
    })

    test_that ("A2: the first one precisely specified", {
        from    <- allOnsets[1]
        till    <- allOnsets[1] + 0.0001   # will cause an error if equal to total recording time
        aSignal <- readEdfSignals (hdr, signals=asn, from=from, till=till, fragment=FALSE)
        aFSignal<- readEdfSignals (hdr, signals=asn, from=from, till=till, fragment=TRUE )
        expect_equal(aSignal, aFSignal)
        expect_equal (aSignal$from , from)
        expect_equal (aSignal$till , till)
        annotsInRange <- (from <= allOnsets) & (allOnsets < till)
        expect_equal (aSignal$annotations, allAnnots[annotsInRange, ])
    })

    test_that ("A3: an empty range precisely specified", {
        if (nOnsets == 1) {
            if (allOnsets[1]) {
                from <- 0
                till <- allOnsets[1] * 0.9999
            } else {  # only one annotation at onset 0
                from <- 1e-8
                till <- Inf
            }
        } else {
            n   <- sample.int(nOnsets-1, 1)
            o1  <- allOnsets[n]             # onset should not before record start
            o2  <- allOnsets[n+1]
            if (o2 >= totalPeriod) o2 <- totalPeriod    # onset may be after totalPeriod
            delta <- abs(o2 - o1) * 0.0001
            from <- o1 + delta
            till <- o2 - delta
        }
        annotsInRange <- (from <= allOnsets) & (allOnsets < till)
        aSignal <- readEdfSignals(hdr, signals=asn, from=from, till=till, fragment=FALSE)
        aFSignal<- readEdfSignals(hdr, signals=asn, from=from, till=till, fragment=TRUE )
        expect_equal(aSignal, aFSignal)
        expect_equal (aSignal$from , from)
        expect_equal (aSignal$till , till)
        expect_equal (sum (annotsInRange), 0)
        expect_equal (nrow (aSignal$annotations), 0)
        if (sum (annotsInRange)) {
            expect_equal (aSignal$annotations, allAnnots[annotsInRange, ])
        }
    })
    test_that ("A4: A random range 1", {
        aRange  <- runif(1) * hdr$recordedPeriod
        from    <- runif(1) * (hdr$recordedPeriod - aRange)
        till    <- from + aRange
        aSignal <- readEdfSignals(hdr, signals=asn, from=from, till=till)
        expect_equal (aSignal$from , from)
        expect_equal (aSignal$till , till)
        annotsInRange <- (from <= allOnsets) & (allOnsets <= till)
        expect_equal (aSignal$annotations, allAnnots[annotsInRange, ])
    })
    test_that ("A5: A random range 2", {
        aRange  <- runif(1) * hdr$recordedPeriod
        from    <- runif(1) * (hdr$recordedPeriod - aRange)
        till    <- from + aRange
        aSignal <- readEdfSignals(hdr, signals=asn, from=from, till=till, fragment=FALSE)
        aFSignal<- readEdfSignals(hdr, signals=asn, from=from, till=till, fragment=TRUE )
        expect_equal(aSignal, aFSignal)
        expect_equal (aSignal$from , from)
        expect_equal (aSignal$till , till)
        annotsInRange <- (from <= allOnsets) & (allOnsets <= till)
        expect_equal (aSignal$annotations, allAnnots[annotsInRange, ])
    })
    test_that ("A6: A random range 3", {
        aRange  <- runif(1) * hdr$recordedPeriod
        from    <- runif(1) * (hdr$recordedPeriod - aRange)
        till    <- from + aRange
        aSignal <- readEdfSignals(hdr, signals=asn, from=from, till=till, fragment=FALSE)
        aFSignal<- readEdfSignals(hdr, signals=asn, from=from, till=till, fragment=TRUE )
        expect_equal(aSignal, aFSignal)
        expect_equal (aSignal$from , from)
        expect_equal (aSignal$till , till)
        annotsInRange <- (from <= allOnsets) & (allOnsets <= till)
        expect_equal (aSignal$annotations, allAnnots[annotsInRange, ])
    })
}

testOSignalsFile <- function (fileNo) {
    hdr <- sHdrs[[fileNo]]
    cat (sFns[fileNo], ": testOSignalsFile: isContinuous=", hdr$isContinuous ,'\n')

    isOSignal   <- !hdr$sHeaders$isAnnotation
    nOSignals   <- sum (isOSignal)
    if (!nOSignals) return

    # select an ordinary signal
    osn <- sample.int(nOSignals, 1)
    osn <- which(isOSignal)[osn]
    signalOsn   <-readEdfSignals (hdr, signals=osn, fragment=FALSE)
    wholeS      <- signalOsn$signal
    nRecords    <- hdr$nRecords
    recordD     <- hdr$recordDuration
    sPerRec     <- hdr$sHeaders$samplesPerRecord[osn]
    sLength     <- length(signalOsn$signal)
    sRate       <- hdr$sHeaders$sRate[osn]
    totalPeriod <- signalOsn$totalPeriod
    sInterval   <- 1/sRate

    sStartRT    <- signalOsn$start
    sStartRS    <- signalOsn$fromSample

    #   The signal samples are counted from 1.
    #   However, as time and duraration are expressed as doubles, time starts at t = 0.
    #   Consequently, with sample rate sRate
    #          the number of samples in [0,t) = ceiling (t * sRate),
    #          the number of samples in [0,t] = floor   (t * sRate) + 1, and
    #          the recording period = [0, s1/sRate) or [0, (s1 - 1) / sRate].
    #   the relation between the edf header paremeter recordDuration and
    #   the number of samples in a record is therefore:
    #          samples in [0:recordDuration) = samples per record
    #   For consistency, the from, till parameters are also interpreted as [from, till)
    #   So fromS, the first sample to read, will be ceiling (from * sRate) + 1

    # at the very beginning
    that    <- "B1: First sample [0:.5)"
    from    <- sStartRT
    till    <- from + sInterval /2
    fromS   <- sStartRS
    part    <- wholeS[1]
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }

    that    <- "B2: First sample [0:1)"
    from    <- sStartRT
    till    <- from + sInterval
    fromS   <- sStartRS
    part    <- wholeS[1]
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }

    that    <- "B3: First 2 samples [0:1.0001)"
    from    <- sStartRT
    till    <- from + sInterval *1.0001
    fromS   <- sStartRS
    part    <- wholeS[1:2]
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }

    that    <- "B4: No first sample [0:0)"
    from    <- sStartRT
    till    <- sStartRT
    fromS   <- NA
    part    <- NULL
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }

    # at the first record boundary
    that    <- "F1: End of record 1: no last sample [end-.9999:end)"
    from    <- sStartRT + recordD - sInterval*.9999
    till    <- sStartRT + recordD
    fromS   <- NA
    part    <- NULL
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }

    that    <- "F2: End of record 1: last sample [end-1:end)"
    from    <- sStartRT + recordD - sInterval
    till    <- sStartRT + recordD
    fromS   <- sStartRS + sPerRec - 1
    part    <- wholeS[sPerRec:sPerRec]
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }

    that    <- "F3: End of record 1: last sample + next [end-1.001:end+0.0001)"
    from    <- sStartRT + recordD - sInterval*1.0001
    till    <- sStartRT + recordD + sInterval*0.0001
    fromS   <- sStartRS + sPerRec - 1
    part    <- wholeS[sPerRec:(sPerRec+1)]
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }

    that    <- "F4: End of 1 record: 1-st sample rec 2 [end:end+1)"
    from    <- sStartRT + recordD
    till    <- sStartRT + recordD + sInterval
    fromS   <- sStartRS + sPerRec
    part    <- wholeS[(sPerRec+1):(sPerRec+1)]
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }

    that    <- "F5: no last sample [end:end)"
    from    <- sStartRT + recordD
    till    <- sStartRT + recordD
    fromS   <- NA
    part    <- NULL
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }


    # the last record
    that    <- "E1: End of recording: no last sample [end-.999:end)"
    from    <- totalPeriod - sInterval*0.9999
    till    <- Inf
    fromS   <- NA
    part    <- NULL
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }

    that    <- "E2: End of recording: last sample [end-1:end)"
    from    <- totalPeriod - sInterval
    till    <- totalPeriod
    fromS   <- sStartRS + sLength - 1
    part    <- wholeS[sLength:sLength]
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }

    that    <- "E3: End of recording: last 2 samples [end-1:end+100)"
    from    <- totalPeriod - sInterval*2
    till    <- totalPeriod + sInterval*100
    fromS   <- sStartRS + sLength - 2
    part    <- wholeS[(sLength-1):sLength]
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }

    that    <- "E4: End of recording: last 36 samples [end-36:inf)"
    from    <- totalPeriod - 36 * sInterval
    till    <- Inf
    fromS   <- sStartRS + sLength - 36
    part    <- wholeS[(sLength-35):sLength]
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }

    # some random ranges
    that    <- "R1: An 1 recordDuration cPart with random start"
    maxFrom <- recordD * (nRecords -1)
    from1   <- runif(1)*maxFrom
    fromS1  <- ceiling (sRate * from1) +1
    from    <- sStartRT + from1
    till    <- from + recordD
    fromS   <- sStartRS + fromS1 -1
    part <- wholeS[fromS1:(fromS1+sPerRec-1)]
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }

    that    <- "R2: An 3 recordDurations cPart with random start"
    maxFrom <- recordD * (nRecords - 3)
    from1   <- runif(1)*maxFrom
    fromS1  <- ceiling (sRate * from1) +1
    from    <- sStartRT + from1
    till    <- from + 3 * recordD
    fromS   <- sStartRS + fromS1 -1
    part    <- wholeS[fromS1:(fromS1+sPerRec*3-1)]
    testCPart (that, hdr, osn, from, till, fromS, part)
    if (!hdr$isContinuous) {
        testFPart (that, hdr, osn, from, till, fromS, part)
    }
}

testCPart <- function (that, hdr, osn, from, till, fromS, wPart) {
    that <- paste ("testCPart: ", that, sep='')
    # cat (that, '\n')
    cPart <- readEdfSignals(hdr, signals=osn, from=from, till=till, fragments = FALSE)
    test_that (that, {
        expect_equal (cPart$from , from)
        expect_equal (cPart$till , till)
        if (is.null(wPart)) {
            expect_equal (length(cPart$signal), 0)
            expect_equal (cPart$fromSample, as.numeric(NA))
            expect_equal (fromS, NA)
        } else {
            expect_equal (cPart$fromSample, fromS)
            expect_equal (cPart$signal, wPart)
        }
    })
}

testFPart <- function (that, hdr, osn, from, till, fromS, wPart) {
    that <- paste ("testFPart: ", that, sep='')
    # cat (that, '\n')
    fSig <- readEdfSignals(hdr, signals=osn, from=from, till=till, fragments=TRUE)
    test_that (that, {
        expect_equal (fSig$till, till)
        expect_equal (fSig$from, from)
    })
    fragments <- fSig$fragments
    if (is.null(wPart) || sum (!is.na(wPart))==0) {
        test_that (that, {
            expect_equal (length(fragments), 0)    # i.e. no fragments
        })
    } else {
        test_that (that, {
            expect_gt (length (fragments), 0)
        })
        wp <- wPart[!is.na(wPart)]
        wPFromFS <- 1
        for (fn in 1:length(fragments)) {
            frag    <- fragments[[fn]]
            # fromRS  <- frag$firstSample
            fSignL  <- length (frag$signal)
            test_that (that, {
                expect_equal(frag$signal, wp[1:fSignL])
            })
            wp <- wp[(fSignL+1):length(wp)]
        }
    }
}

set.seed (20160101)

testAll <- function () {
    testASignalsFile (1)
    testASignalsFile (2)
    testASignalsFile (3)
    testASignalsFile (4)
    testOSignalsFile (1)
    testOSignalsFile (2)
    testOSignalsFile (3)
    testOSignalsFile (4)
}

testAll()
