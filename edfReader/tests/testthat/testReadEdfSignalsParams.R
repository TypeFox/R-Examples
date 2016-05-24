#
# Purpose   :    Test the use of readEdfSignals parameters
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
#   Mar16 - Revised (somewahat more generic & support added for :
#           files with more then one annotation signal and first record onset != 0
# ------------------------------------------------------------------------------
require (testthat)
require (edfReader)

context ("Compare reading whole files with the export from EDFBrowser.")

libDir <- paste (system.file("extdata", package="edfReader"), '/', sep='')

sFns <- c('edfPlusC.edf',
          'bdfPlusC.bdf',          # same signals as in edfPlusC.edf, but 24 bits
          'edfPlusD.edf',          # wextracted from 'edfPlusC.edf'
          'edfAnnonC.edf')         # with 2 annotation signals, to test annotations only'

# edfPlusC.edf and bdfPlusC.bdf are truncated versions of the test_generator_2 test file
# the original files can be found at: http://www.teuniz.net/edf_bdf_testfiles
# edfAnnonC.edf  is a truncated version of the test_generator8.edf test file
# test_generator8.edf has been received from Teunis van Beelen via private communications

sFFns <- character (length = length(sFns))
for (i in 1:length(sFns)) sFFns[i] <- paste (libDir, sFns[i], sep='')

sHdrs <- vector (mode='list', length = length(sFns))
for (i in 1:length(sFns)) sHdrs[[i]] <- readEdfHeader(sFFns[i])

#                            Test ranges
# ------------------------------------------------------------------------------
testMinMax <- function  (sFn) {
    cat (sFns[sFn], ': min/max tests\n')
    hdr         <- sHdrs[[sFn]]

    #  digital range
    that <- "the digital/physical Min value must be smaller than the digital Max value"

    hdr <- sHdrs[[1]]
    hdr$sHeaders$digitalMin[2]  <-  hdr$sHeaders$digitalMax[2]
    expect_error (readEdfSignals (hdr), "Illegal digital min/max, use physical=FALSE",
                  info="the digitalMin should be less than digitalMax")

    #  physical range
    hdr <- sHdrs[[1]]
    hdr$sHeaders$physicalMin[3]  <-  hdr$sHeaders$physicalMax[3]
    expect_error (readEdfSignals (hdr), "Illegal physical min/max, use physical=FALSE",
                  info="the physicalMin should be less than physicalMax")

    # both ranges
    hdr <- sHdrs[[1]]
    hdr$sHeaders$digitalMin[2]   <-  hdr$sHeaders$digitalMax[2]
    hdr$sHeaders$physicalMin[3]  <-  hdr$sHeaders$physicalMax[3]
    expect_error (readEdfSignals (hdr), "Illegal digital/physical min/max, use physical=FALSE",
                  info="both the digital/physicalMin should be less than digital/physicalMax")
}

#                           Test label selections
# ------------------------------------------------------------------------------
testSignalSelection <- function (sFn) {
    cat (sFns[sFn], ': signal selection tests\n')
    that <- 'signal label selections'
    hdr         <- sHdrs[[sFn]]
    nHSignals   <-  hdr$nSignals
    sAll <- readEdfSignals(hdr,simplify = FALSE)
    allSignalsN <- length (sAll)

    s <- readEdfSignals(hdr, signals='Annotations', simplify = FALSE)
    SASignalsL  <- sapply (s,  function(X) {X$isAnnotation})
    expect_equal (sum(!SASignalsL), 0,
                  info="the only signals read should be annotations signals")

    s <- readEdfSignals(hdr, signals='Ordinary', simplify = FALSE)
    SOSignalsL  <- sapply (s,  function(X) {!X$isAnnotation})
    expect_equal (sum(!SOSignalsL), 0,
                 info="the only signals read should be oridinary signals")
    expect_equal(sum(SASignalsL)+sum(SOSignalsL), allSignalsN,
                info="there shouldn't be an annotation signal")

    s <- readEdfSignals(hdr, signals=hdr$nSignals)
    expect_equal (s$signalNumber, hdr$nSignals)
    oob <-  hdr$nSignals + 1
    expect_error(readEdfSignals(hdr, signals=oob), paste ("Signal number out of bound:", oob))
    expect_error(readEdfSignals(hdr, signals=c(100, 200)), "Signal numbers out of bound: 100 200")

    expect_error (readEdfSignals(hdr, signals='Bla'), "Unkown signal designation Bla")
    expect_error(readEdfSignals(hdr, signals=""), "Unkown signal designation")

    if (hdr$fileType == 'BDF') {
        aLabel <- "BDF Annotations"
        aName1 <- "BDF Annotations-1"
    } else {
        aLabel <- "EDF Annotations"
        aName1 <- "EDF Annotations-1"
    }

    s <- readEdfSignals(hdr)
    if (aName1 %in% names(s)) {
        asn1 <- which.max (names(s)==aName1)
        s <- readEdfSignals(hdr, signals=c(asn1, as.character(asn1), aName1, aLabel), simplify=FALSE)
        expect_equal (length(s), 1)
        expect_equal (s[[asn1]]$signalNumber,  which(hdr$sHeaders$isAnnotation))
    }

    s <- readEdfSignals(hdr)
    oSL     <- sapply (s,  function(X) {!X$isAnnotation})
    oS1n    <- which.max(oSL)
    oS1c    <- as.character(oS1n)
    oName1  <- names(s)[oS1n]
    oLabel1 <- s[[oS1n]]$label
    oSall   <- readEdfSignals(hdr, signals =c(oS1n, oS1c, oName1, oLabel1), simplify=FALSE)
    oS1n    <- readEdfSignals(hdr, signals= oS1n,    simplify=FALSE)
    oS1c    <- readEdfSignals(hdr, signals= oS1c,    simplify=FALSE)
    oSl1    <- readEdfSignals(hdr, signals= oLabel1, simplify=FALSE)
    oSn1    <- readEdfSignals(hdr, signals= oName1,  simplify=FALSE)
    oSn1S   <- readEdfSignals(hdr, signals= oName1)
    expect_equal (oS1n[[1]], oSall[[1]])
    expect_equal (oS1n[[1]], oS1c [[1]])
    expect_equal (oS1n[[1]], oSn1 [[1]])
    expect_equal (oS1n[[1]], oSl1 [[1]])
    expect_equal (oS1n[[1]], oSn1S)
    if (oName1 == oLabel1) {
        expect_equal (length(oSl1),  1)
        expect_equal (length(oSall), 1)
    } else {
        expect_gt (length(oSl1),  1)
        expect_gt (length(oSall), 1)
    }
}

testAllFiles <- function () {
    testMinMax(1)
    testMinMax(2)
    testMinMax(3)
    testMinMax(4)
    testSignalSelection (1)
    testSignalSelection (2)
    testSignalSelection (3)
    testSignalSelection (4)
}

testAllFiles()
