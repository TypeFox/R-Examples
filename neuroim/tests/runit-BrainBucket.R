# TODO: Add comment
# 
# Author: brad
###############################################################################


test.BrainBucketSource <- function() {
	bsource <- BrainBucketSource("data/statout+orig.HEAD")
	checkTrue(!is.null(bsource))
	checkEquals(length(bsource@metaInfo@label), 25)
}

test.loadBucket <- function() {
	buck <- loadBucket("data/statout+orig.HEAD")
}

test.eachVolume <- function() {
	buck <- loadBucket("data/statout+orig.HEAD")
	mean.vol1 <- eachVolume(buck, mean, na.rm=TRUE)
	checkEquals(length(mean.vol1), 25)
}
