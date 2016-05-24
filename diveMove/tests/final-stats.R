library(diveMove)

zz <- system.file(file.path("data", "dives.csv"),
                  package="diveMove", mustWork=TRUE)
(sealX <- readTDR(zz, concurrentCols=4:6, sep=";",
                  na.strings="", as.is=TRUE))
(dcalib <- calibrateDepth(sealX, zoc.method="offset", offset=3))
tdrstats <- diveStats(dcalib)
head(tdrstats)
head(stamps <- stampDive(dcalib))
(att <- timeBudget(dcalib, FALSE))
(att <- timeBudget(dcalib, TRUE))

(sealX <- readTDR(zz, speed=TRUE, concurrentCols=4:6,
                  sep=";", na.strings="", as.is=TRUE))
(dcalib <- calibrateDepth(sealX, zoc.method="offset", offset=3))
(vcalib <- calibrateSpeed(dcalib, z=1))
tdrstats <- diveStats(vcalib)
head(tdrstats)
head(stamps <- stampDive(vcalib))
(att <- timeBudget(vcalib, FALSE))
(att <- timeBudget(vcalib, TRUE))


###_ + Emacs local variables
## Local variables:
## allout-layout: (+ : 0)
## End:
