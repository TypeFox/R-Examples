
d <- tracks( matrix(c(0, 8,  
		10, 9, 
		20, 7, 
		30, 7,
		40, 6, 
		50, 5), ncol=2, byrow=T ) )

test_that("all track measures work", {
	expect_equal(trackLength(d[[1]]), 5)
	expect_equal(duration(d[[1]]), 50)
	expect_equal(speed(d[[1]]), 1/10)
	expect_equal(displacement(d[[1]]), 3)
	expect_equal(squareDisplacement(d[[1]]), 9)
	expect_equal(maxDisplacement(d[[1]]), 3)
	expect_equal(displacementRatio(d[[1]]), 1)
	expect_equal(outreachRatio(d[[1]]), 3/5)
	expect_equal(overallAngle(d[[1]]), pi)
	expect_equal(meanTurningAngle(d[[1]]), pi/2)
	expect_equal(overallDot(d[[1]]), -1)
	expect_equal(asphericity(d[[1]]), NaN)
	expect_equal(hurstExponent(d[[1]]), 0)
	expect_equal(fractalDimension(d[[1]]), 1.415037499278843480255)
})


test_that("aggregation of track measures works", {
	expect_equal(aggregate(d,trackLength)[,2], c(1,2,8/3,4,5))
	expect_equal(aggregate(d,duration)[,2], 10*c(1,2,3,4,5))
	expect_equal(aggregate(d,duration)[,2], 10*c(1,2,3,4,5))
})