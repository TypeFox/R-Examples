# TODO: Add comment
# 
# Author: brad
###############################################################################


test.DenseBrainVector <- function() {
	dat <- array(0, c(64,64,64,4))
	spc <- BrainSpace(c(64,64,64,4))
	bv <- DenseBrainVector(dat, spc)
	checkTrue(!is.null(bv))
	checkEquals(bv[1,1,1,1], 0)
	checkEquals(bv[64,64,64,4], 0)
	checkEquals(dim(bv), c(64,64,64,4))
}

test.DenseBrainVector.concat <- function() {
	dat <- array(0, c(64,64,64,4))
	spc <- BrainSpace(c(64,64,64,4))
	bv1 <- DenseBrainVector(dat, spc)
	bv2 <- DenseBrainVector(dat, spc)
	
	bv3 <- concat(bv1, bv2)
	checkTrue(inherits(bv3, "BrainVector"))
	checkEquals(dim(bv3), c(64,64,64,8))
	
	bv4 <- concat(bv1,bv2, bv1, bv3)
	checkTrue(inherits(bv4, "BrainVector"))
	checkEquals(dim(bv4), c(64,64,64,20))
	checkEquals(bv4[1,1,1,1],0)
	
}

test.DenseBrainVector.takeVolume <- function() {
	dat <- array(0, c(64,64,64,4))
	spc <- BrainSpace(c(64,64,64,4))
	bv1 <- DenseBrainVector(dat, spc)
	bv2 <- DenseBrainVector(dat, spc)
	
	bv3 <- concat(bv1, bv2)
	
	vol1 <- takeVolume(bv3, 1)
	checkEquals(dim(vol1), c(64,64,64))
	
	vec1 <- takeVolume(bv3, 1:2)
	checkTrue(inherits(vec1, "list"))
	checkEquals(dim(vec1[[1]]), c(64,64,64))
	
	vec2 <- takeVolume(bv3, 1:2, merge=TRUE)
	checkTrue(inherits(vec2, "BrainVector"))
	checkEquals(dim(vec2), c(64,64,64,2))
	
}





test.BrainVector.eachVolume <- function() {
	dat <- array(rnorm(64*64*64*4), c(64,64,64,4))
	spc <- BrainSpace(c(64,64,64,4))
	bv1 <- DenseBrainVector(dat, spc)
	
	mean.vol1 <- eachVolume(bv1, mean)
	checkEquals(length(mean.vol1), 4)
	
	mean.vol2 <- eachVolume(bv1, withIndex=TRUE, FUN=function(slice, i) mean(slice))
	checkEquals(length(mean.vol2), 4)
	checkEquals(mean.vol1, mean.vol2)
	
	checkEquals(unlist(mean.vol1), apply(bv1, 4, mean))
}


test.BrainVector.series <- function() {
	dat <- array(rnorm(64*64*64*4), c(64,64,64,4))
	spc <- BrainSpace(c(64,64,64,4))
	bv1 <- DenseBrainVector(dat, spc)
	checkEquals(series(bv1, 1,1,1), bv1[1,1,1,])
	
    mat <- rbind(c(1,1,1), c(2,2,2), c(3,3,3))
	
	r1 <- apply(mat, 1, function(i) { series(bv1, i[1], i[2], i[3]) })
			
	r2 <- series(bv1, mat)
	
	checkEquals(r1, r2)
	
	
}


test.BrainVector.as.matrix <- function() {
	dat <- array(rnorm(64*64*64*4), c(64,64,64,4))
	spc <- BrainSpace(c(64,64,64,4))
	bv1 <- DenseBrainVector(dat, spc)
	mat <- as.matrix(bv1)
	
	ind <- 1:(64*64*64)
	mat2 <- t(do.call(cbind, lapply(ind, function(i) series(bv1, i))))
	
	checkEquals(mat, mat2)
	
}

test.BrainVector.roundtrip.io <- function() {
	bvec <- loadVector("data/qrscan01.nii.gz", indices=1:4)
	template <- takeVolume(bvec,1)
	
	mask.idx <- sort(sample(1:length(template), 1000))
	vals <- rnorm(length(mask.idx))
	bv <- BrainVolume(vals, space(template), indices=mask.idx)
	fname <- paste(tempfile(), ".nii", sep="")
	writeVolume(bv,fname)
	bv2 <- loadVolume(fname)
	
	checkEquals(dim(bv2), dim(bv))
	checkEqualsNumeric(trans(bv2), trans(bv), tol=.0001)
	
}

test.SparseBrainVector.series <- function() {
	
	spvec <- loadVector("data/qrscan01.nii.gz", indices=1:4, mask=rep(TRUE, 96*96*26))
	bvec <- loadVector("data/qrscan01.nii.gz", indices=1:4)
	
	voxmat <- rbind(c(10,10,10), c(20,20,10), c(30,30, 10), c(40,40,10), c(50,50,10))
	
	checkEquals(series(bvec, voxmat), series(spvec, voxmat))
}
	

test.SparseBrainVector.roundtrip.io <- function() {
	bvec <- loadVector("data/qrscan01.nii.gz", indices=1:4, mask=rep(TRUE, 96*96*26))
	template <- takeVolume(bvec,1)
	
	mask.idx <- sort(sample(1:length(template), 1000))
	vals <- rnorm(length(mask.idx))
	bv <- BrainVolume(vals, space(template), indices=mask.idx)
	fname <- paste(tempfile(), ".nii", sep="")
	writeVolume(bv,fname)
	bv2 <- loadVolume(fname)
	
	checkEquals(dim(bv2), dim(bv))
	checkEqualsNumeric(trans(bv2), trans(bv), tol=.0001)
	
}

test.SparseBrainVector.as.sparse <- function() {
	dat <- array(rnorm(64*64*64*4), c(64,64,64,4))
	spc <- BrainSpace(c(64,64,64,4))
	bv1 <- DenseBrainVector(dat, spc)
	svec <- as.sparse(bv1, c(1,100,1000))
	checkEquals(dim(svec), c(64,64,64,4))
	
	
}

test.SparseBrainVector <- function() {
	dat <- array(rnorm(64*64*64*4), c(64,64,64,4))
	spc <- BrainSpace(c(64,64,64,4))
	tmp <- rnorm(64*64*64)
	mask <- tmp < 1000
	mask <- LogicalBrainVolume(mask, dropDim(spc))
	bvec <- SparseBrainVector(dat, spc, mask)
	
	checkEquals(dim(bvec), dim(dat))
	checkEquals(dat[1,1,1,1], bvec[1,1,1,1])
	checkEquals(dat[1,1,1,1], bvec[1,1,1,1])
	
}

test.SparseBrainVector.concat <- function() {
	dat <- array(0, c(64,64,64,4))
	spc <- BrainSpace(c(64,64,64,4))
	tmp <- rnorm(64*64*64)
	mask <- tmp > .8
	mask <- LogicalBrainVolume(mask, dropDim(spc))
	
	bv1 <- SparseBrainVector(dat, spc, mask)
	
	bv2 <- concat(bv1, bv1)
	checkTrue(inherits(bv2, "BrainVector"))
	checkEquals(dim(bv2), c(64,64,64,8))
	
	bv3 <- concat(bv1,bv2, bv1, bv2)
	checkTrue(inherits(bv3, "BrainVector"))
	checkEquals(dim(bv3), c(64,64,64,24))
	#checkEquals(bv4[1,1,1,1],0)
}
	


	


