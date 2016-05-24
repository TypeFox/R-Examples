### File R/sisalData.R
### This file is part of the sisal package for R.
###
### Copyright (C) 2015 Aalto University
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### A copy of the GNU General Public License is available at
### http://www.r-project.org/Licenses/

.MD5laser        <- "2a4071c2b9bed90a09c65f3f79c1f692"
.MD5laser.cont   <- "4cf1b2ccd6a0e3d16cb5519eeaca0b8e"
.MD5poland.learn <- "707caa3c7711ff7091b397c36d361c1e"
.MD5poland.test  <- "ce13f9204de00a362e53c1431370b229"

sisalData <- function(dataset = c("laser", "laser.cont", "poland"),
                      verify = TRUE) {

    dataset2 <- match.arg(dataset)
    stopifnot(is.logical(verify), !is.na(verify), length(verify) == 1)

    if (dataset2 == "laser") {
        foo <- read.table("http://www-psych.stanford.edu/~andreas/Time-Series/SantaFe/A.dat")[[1]]
        if (verify && !identical(digestVector(foo), .MD5laser)) {
            stop("remote file has changed")
        }
        foo
    } else if (dataset2 == "laser.cont") {
        foo <- read.table("http://www-psych.stanford.edu/~andreas/Time-Series/SantaFe/A.cont")[[1]]
        if (verify && !identical(digestVector(foo), .MD5laser.cont)) {
            stop("remote file has changed")
        }
        foo
    } else if (dataset2 == "poland") {
        fname <- tempfile()
        if (download.file(url="http://research.ics.aalto.fi/eiml/datasets/PolandElectricity.zip",
                          destfile=fname, mode="wb", quiet=TRUE) != 0) {
            stop("failed to download file")
        }
        on.exit(unlink(fname))
        dname <- tempfile()
        if (!dir.create(dname)) {
            stop("failed to create directory")
        }
        on.exit(unlink(dname, recursive=TRUE), add=TRUE)
        learnFile <- "ElectricFixed.mat"
        testFile <- "ElectricFixedTest.mat"
        if (length(unzip(zipfile=fname, exdir=dname,
                         files=c(learnFile, testFile))) != 2) {
            stop()
        }
        poland.learn <-
            drop(suppressMessages(suppressWarnings(readMat(file.path(dname,
                                                                     learnFile))[[1]])))
        if (verify &&
            !identical(digestVector(poland.learn), .MD5poland.learn)) {
            stop("remote file has changed")
        }
        poland.test <-
            drop(suppressMessages(suppressWarnings(readMat(file.path(dname,
                                                                     testFile))[[1]])))
        if (verify &&
            !identical(digestVector(poland.test), .MD5poland.test)) {
            stop("remote file has changed")
        }
        list(learn=poland.learn, test=poland.test)
    }
}
