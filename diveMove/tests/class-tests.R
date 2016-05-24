library(diveMove)

###_ + Reading and as.data.frame ------------------------------------------
zz <- system.file(file.path("data", "dives.csv"),
                  package="diveMove", mustWork=TRUE)
sealX <- readTDR(zz, speed=TRUE, sep=";", na.strings="", as.is=TRUE)
stopifnot(is(sealX, "TDR"), is(sealX, "TDRspeed"))

sealX <- readTDR(zz, sep=";", na.strings="", as.is=TRUE)
stopifnot(is(sealX, "TDR"), !is(sealX, "TDRspeed")) #speed=FALSE default

sealDat <- as.data.frame(sealX)
stopifnot(is.data.frame(sealDat))
sealX <- readTDR(zz, concurrentCols=NULL, sep=";",
                 na.strings="", as.is=TRUE)
## concurrentData should be a 0-column data.frame
stopifnot(ncol(getCCData(sealX)) == 0)

sealX <- readTDR(zz, subsamp=10, speed=TRUE, concurrentCols=6,
                 sep=";", na.strings="", as.is=TRUE)
stopifnot(is(sealX, "TDR"), is(sealX, "TDRspeed"))
stopifnot(getDtime(sealX) == 10)
sealX <- readTDR(zz, subsamp=10, speed=TRUE, concurrentCols=5:6,
                 sep=";", na.strings="", as.is=TRUE)
stopifnot(ncol(getCCData(sealX)) == 2)

sealDat <- as.data.frame(sealX)
sealX <- createTDR(time=sealDat$time, depth=sealDat$depth,
                   concurrentData=sealDat[, 3:ncol(sealDat)],
                   dtime=sealX@dtime, file=sealX@file)
sealX.speed <- createTDR(time=sealDat$time, depth=sealDat$depth,
                         concurrentData=sealDat[, 3:ncol(sealDat)],
                         speed=TRUE, dtime=sealX@dtime, file=sealX@file)
stopifnot(is(sealX, "TDR"), is(sealX.speed, "TDRspeed"))

###_ + Accessors ----------------------------------------------------------
head(tt <- getTime(sealX))
head(dd <- getDepth(sealX))
head(ss <- getSpeed(sealX.speed))
head(cc <- getCCData(sealX))
head(cc <- getCCData(sealX, "speed"))
stopifnot(identical(length(tt), length(dd), length(ss), ncol(cc)))

getFileName(sealX)
getDtime(sealX)

###_ + Replacements -------------------------------------------------------
sll <- length(getSpeed(sealX.speed))
speed.new <- rnorm(sll)
speed(sealX.speed) <- speed.new
stopifnot(identical(speed.new, getSpeed(sealX.speed)))
head(getSpeed(sealX.speed))

depth.new <- rnorm(getDepth(sealX))
depth(sealX) <- depth.new
stopifnot(identical(depth.new, getDepth(sealX)))
sealX <- createTDR(time=sealDat$time, depth=sealDat$depth,
                   concurrentData=sealDat[, 3:ncol(sealDat)],
                   dtime=sealX@dtime, file=sealX@file)
depth.new <- rnorm(length(getDepth(sealX)))
depth(sealX) <- depth.new
stopifnot(identical(depth.new, getDepth(sealX)))
head(getDepth(sealX))


###_ + Emacs local variables
## Local variables:
## allout-layout: (+ : 0)
## End:
