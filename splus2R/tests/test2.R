# This file is for checking expectWarnings() and expectStop()

library("splus2R")

testFile <- "warnings.t"
do.test(testFile)
do.test(testFile, local=TRUE)

testFile <- "stops.t"
do.test(testFile)
do.test(testFile, verbose=TRUE)
