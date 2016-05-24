# TODO: Add comment
# 
# Author: brad
###############################################################################


test.makeKernel <- function() {
	mask <- loadVolume("data/global_mask.nii")
	kern <- makeKernel(c(5,5,3), spacing(mask), function(v) dnorm(v, sd=3))
	checkTrue(!is.null(kern))
	checkEquals(dim(kern$coordmat), c(75,3))
}

test.indexToGrid <- function() {
	mask <- loadVolume("test_data/global_mask.nii")
	idx <- which(mask>0)
	grid <- indexToGrid(mask, idx)
	checkEquals(nrow(grid), length(idx))
}


test.gridToIndex <- function() {
	mask <- loadVolume("test_data/rscan001_mask.nii.gz")
	bvec <- loadVector("test_data/rscan001.nii.gz", mask=mask)
	
	MAXIND <- prod(dim(mask))
	
	idx <- which(mask > 0)
	gridmat <- indexToGrid(mask, idx)   
	
	kernel <- makeKernel(c(5,5,3), spacing(mask), function(v) dnorm(v, sd=3))
	kwts <- kernel$wts
	indmat <- as.matrix(kernel$indmat)
	
}



