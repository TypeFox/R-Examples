context("CCA Single")
test_that("cca single works for full matrix", {
	expect_that(cca.single(matrix(1:12,nrow=3), s=1, x=3,y=1,mode=1), equals(matrix(rep(1,12),nrow=3)))
	expect_that(cca.single(matrix(1:12,nrow=3), s=1, x=3,y=1,mode=2), equals(matrix(rep(1,12),nrow=3)))
	expect_that(cca.single(matrix(1:12,nrow=3), s=1, x=3,y=1,mode=3), equals(matrix(rep(1,12),nrow=3)))
		})
test_that("cca single works for half matrix", {
	half.matrix <- matrix(rep(c(1,0),each=6),nrow=3)
	expect_that(cca.single(half.matrix, s=1, x=1,y=3,mode=1), equals(matrix(rep(0,12),nrow=3)))
	expect_that(cca.single(half.matrix, s=1, x=1,y=2,mode=1), equals(half.matrix))
		})

test_that("cca single works differently for nearest neighbors (mode=1), shell (mode=2) and radius (mode=3) ", {
	complex.matrix <- matrix(rep(0,60),nrow=12)
	complex.matrix[1,1] <- 1
	corner.only <- complex.matrix
	complex.matrix[2,2] <- 1 #not reached by nearest neighbor, but by shell of size 1
	two.areas <- complex.matrix
	complex.matrix[4,4] <- 1 #not reached by radius 2 neighborhood, but by shell of size 2
	expect_that(cca.single(complex.matrix, s=1, x=1,y=1,mode=1), equals(corner.only))
	expect_that(cca.single(complex.matrix, s=1, x=1,y=1,mode=2), equals(two.areas))
	expect_that(cca.single(complex.matrix, s=2, x=1,y=1,mode=3), equals(two.areas))
	expect_that(cca.single(complex.matrix, s=2, x=1,y=1,mode=2), equals(complex.matrix))
		})
test_that("cca single works equivalent for nearest neighbors (mode=1) and radius with s=1 (mode=3) ", {
	complex.matrix <- matrix(rep(0,60),nrow=12)
	complex.matrix[1,1] <- 1
	corner.only <- complex.matrix
	complex.matrix[2,2] <- 1 #not reached by nearest neighbor, but by shell of size 1
	two.areas <- complex.matrix
	complex.matrix[4,4] <- 1 #not reached by radius 2 neighborhood, but by shell of size 2
	expect_that(cca.single(complex.matrix, s=1, x=1,y=1,mode=3), equals(corner.only))
		})

