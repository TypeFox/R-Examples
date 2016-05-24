


test.DenseBrainVolume <- function() {
	dat <- array(0, c(64,64,64))
	spc <- BrainSpace(c(64,64,64))
	bv <- DenseBrainVolume(dat, spc)
	checkTrue(!is.null(bv))
	checkEquals(bv[1,1,1], 0)
	checkEquals(bv[64,64,64], 0)
	checkEquals(dim(bv), c(64,64,64))
}

test.DenseBrainVolume.indices <- function() {
	dat <- rnorm(100)
	spc <- BrainSpace(c(64,64,64))
	
	indices = seq(1,20000, length.out=100)
	bv <- DenseBrainVolume(dat, spc, indices=indices)
	checkTrue(!is.null(bv))
	checkEquals(dim(bv), c(64,64,64))
}

test.DenseBrainVolume.concat <- function() {
	dat <- array(0, c(64,64,64))
	spc <- BrainSpace(c(64,64,64))
	bv1 <- DenseBrainVolume(dat, spc)
	bv2 <- DenseBrainVolume(dat, spc)
	
	bv3 <- concat(bv1, bv2)
	checkTrue(inherits(bv3, "BrainVector"))
	checkEquals(dim(bv3), c(64,64,64,2))
	
	bv4 <- concat(bv1,bv2, bv1, bv2)
	checkTrue(inherits(bv4, "BrainVector"))
	checkEquals(dim(bv4), c(64,64,64,4))
	checkEquals(bv4[1,1,1,1],0)
	
}

test.loadVolume <- function() {
	vol <- loadVolume("data/global_mask.nii")
	checkTrue(!is.null(vol))
	checkEquals(dim(vol), c(64,64,25))
	
	checkException(loadVolume("data/global_mask.nii", index=5), silent=TRUE)
	
	checkEquals(dim(loadVolume("data/epivector.nii", index=3)), c(64,64,38))
	
}

test.loadVolume.gz <- function() {
	vol <- loadVolume("data/fixef_symbol2_diff_GLT#0_Tstat.nii.gz")
	checkTrue(!is.null(vol))
	checkEquals(dim(vol), c(96,96,26))
	
}

test.loadVolume2.gz <- function() {
	vol <- loadVolume("data/blade_zippedvol.nii.gz")
	checkTrue(sum(vol) == 54866)
	
}



test.vol.arith <- function() {
	vol1 <- loadVolume("data/global_mask.nii")
	vol2 <- loadVolume("data/global_mask.nii")
	
	vol.plus <- vol1 + vol2
	vol.check <- DenseBrainVolume(vol1@.Data + vol2@.Data, space(vol1))
	checkEquals(vol.plus, vol.check)
	
	vol.minus <- vol1 - vol2
	vol.check <- DenseBrainVolume(vol1@.Data - vol2@.Data, space(vol1))
	checkEquals(vol.minus, vol.check)
	
	vol.mult <- vol1 * vol2
	vol.check <- DenseBrainVolume(vol1@.Data * vol2@.Data, space(vol1))
	checkEquals(vol.mult, vol.check)
	
}

test.vol.indexToGrid <- function() {
	vol1 <- loadVolume("data/global_mask.nii")
	
	i <- 65
	checkEquals(indexToGrid(vol1, i), matrix(c(1,2,1), nrow=1))
	checkException(indexToGrid(vol1, 64*64*25 +1), silent=TRUE)
	checkException(indexToGrid(vol1, -1), silent=TRUE)
}

test.LogicalBrainVolume <- function() {
	vol1 <- loadVolume("data/global_mask.nii")
	vol2 <- as(vol1, "LogicalBrainVolume")
	checkTrue(!is.null(vol2))
}


test.vol.eachSlice <- function() {
	vol1 <- loadVolume("data/global_mask.nii")
	mean.slice1 <- eachSlice(vol1, mean)
	checkEquals(length(mean.slice1), 25)
	
	mean.slice2 <- eachSlice(vol1, withIndex=TRUE, FUN=function(slice, i) mean(slice))
	checkEquals(length(mean.slice2), 25)
	checkEquals(mean.slice1, mean.slice2)
	
	checkEquals(unlist(mean.slice1), apply(vol1, 3, mean))
}

test.writeVolume <- function() {
	vol1 <- loadVolume("data/global_mask.nii")
	fname <- paste(tempfile(), ".nii", sep="")
	writeVolume(vol1, fname)
	
	
	vol2 <- loadVolume(fname)
	checkTrue(all(vol1 == vol2))
	checkEquals(vol2@source@metaInfo@dataType, vol1@source@metaInfo@dataType)
	
	checkTrue(identical(space(vol1), space(vol2)))
	
	fname <- paste(tempfile(), ".nii", sep="")
	writeVolume(vol1, fname, dataType="DOUBLE")
	vol3 <- loadVolume(fname)
	checkEquals(vol3@source@metaInfo@dataType, "DOUBLE")
	
}



