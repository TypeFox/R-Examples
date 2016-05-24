### R code from vignette source 'h5.Rnw'

###################################################
### code chunk number 1: h5.Rnw:47-51
###################################################
prettyVersion <- packageDescription("h5")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")
require(h5)
require(highlight)


###################################################
### code chunk number 2: h5.Rnw:158-164
###################################################
f <- h5file("test.h5")
f["testgroup/testset"] <- rnorm(100)
testattr <- LETTERS[round(runif(100, max=26))]
h5attr(f["testgroup/testset"], "testattr") <- testattr
f["testgroup/testset"]
h5close(f)


###################################################
### code chunk number 3: h5.Rnw:213-220
###################################################
f <- h5file("test.h5")
f["testmat"] <- matrix(rep(1L, 6), nrow=3)
f["testmat"][c(1, 3), 2] <- rep(2L, 2)
cbind(f["testmat"], matrix(7:9, nrow=3))
f["testmat"][]
h5unlink(f, "testmat")
h5close(f)


###################################################
### code chunk number 4: h5.Rnw:241-250
###################################################
library(zoo)
datevec <- seq(as.POSIXct("2015-12-01"), as.POSIXct("2016-01-01"), by = "secs")
tsdat <- zoo(matrix(rnorm(length(datevec) * 3), ncol=3), order.by=datevec)
f <- h5file("test.h5", "a")
f["testseries", chunksize=c(86400, 1)] <- cbind(index(tsdat), coredata(tsdat))
h5flush(f)
tssub <- zoo(f["testseries"][1:86400, 2], order.by=as.POSIXct(f["testseries"][1:86400, 1], origin="1970-01-01"))
identical(tssub, tsdat[1:86400, 1, drop=FALSE])
h5close(f)


###################################################
### code chunk number 6: h5.Rnw:276-279
###################################################
f <- h5file("ex-matlab.mat", "r")
dates <- as.Date(f["tseries"][1, 1:3] - 719529)
zoo(t(f["tseries"][2:4, 1:3]), order.by=dates)


###################################################
### code chunk number 8: h5.Rnw:302-305
###################################################
f <- h5file("ex-pandas.h5", "r")
dates <- as.Date(f["testset/axis1"][1:3] - 719163, origin="1970-01-01")
zoo(f["testset/block0_values"][1:3, ], order.by=dates)


###################################################
### code chunk number 9: h5.Rnw:327-348
###################################################
library(zoo)
library(microbenchmark)
# Create data files
datevec <- seq(as.POSIXct("2015-12-01"), as.POSIXct("2016-01-01"), by = "secs")
tsdat <- zoo(matrix(rnorm(length(datevec) * 3), ncol=3), order.by=datevec)
f <- h5file("testbm.h5", "a")
f["testseries", chunksize=c(86400, 1)] <- cbind(index(tsdat), coredata(tsdat))
h5close(f)
save(tsdat, file="test.rda")

readRda <- function() {
  load("test.rda")
  tsdat[1:86400, 1:2]
}
readH5 <- function() {
  f <- h5file("test.h5", "r")
  f["testseries"][1:86400, 1:2]
  h5close(f)
}
bm <- microbenchmark(readRda(), readH5(), times = 10L)
summary(bm)$median[1] / summary(bm)$median[2]


###################################################
### code chunk number 10: h5.Rnw:351-352
###################################################
dummy <- file.remove(c("test.h5", "test.rda", "testbm.h5"))


