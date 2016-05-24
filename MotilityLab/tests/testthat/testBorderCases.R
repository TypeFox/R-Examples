
d <- tracks( matrix(c(0, 8,  
		10, 9, 
		20, 7, 
		30, 7,
		40, 6, 
		50, 5), ncol=2, byrow=T ),
	matrix( c(0,0), nrow=1 ) )


test_that("all track measures work on zero-length tracks", {
	expect_equal(trackLength(matrix()), 0)
	expect_equal(duration(matrix()), 0)
	expect_equal(speed(matrix()), NaN)
	expect_equal(displacement(matrix()), 0)
	expect_equal(squareDisplacement(matrix()), 0)
	expect_equal(maxDisplacement(matrix()), 0)
	expect_equal(displacementRatio(matrix()), NaN)
	expect_equal(outreachRatio(matrix()), NaN)
	expect_equal(overallAngle(matrix()), 0)
	expect_equal(meanTurningAngle(matrix()), NaN)
	expect_equal(overallDot(matrix()), NaN)
	expect_equal(asphericity(matrix()), 1)
	expect_equal(hurstExponent(matrix()), NA)
})

test_that("all track measures work on tracks of length one", {
	expect_equal(trackLength(d[[2]]), 0)
	expect_equal(duration(d[[2]]), 0)
	expect_equal(speed(d[[2]]), NaN)
	expect_equal(displacement(d[[2]]), 0)
	expect_equal(squareDisplacement(d[[2]]), 0)
	expect_equal(maxDisplacement(d[[2]]), 0)
	expect_equal(displacementRatio(d[[2]]), NaN)
	expect_equal(outreachRatio(d[[2]]), NaN)
	expect_equal(overallAngle(d[[2]]), 0)
	expect_equal(meanTurningAngle(d[[2]]), NaN)
	expect_equal(overallDot(d[[2]]), NaN)
	expect_equal(asphericity(d[[2]]), 1)
	expect_equal(hurstExponent(d[[2]]), NA)
})


test_that("aggregation of track measures works", {
	expect_equal(aggregate(d,trackLength)[,2], c(1,2,8/3,4,5))
	expect_equal(aggregate(d,duration)[,2], 10*c(1,2,3,4,5))
})