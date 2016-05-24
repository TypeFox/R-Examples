
set.seed(2345)

test_that("Track simulation", {
	expect_equal(nrow(brownianTrack(100,1)), 101)
	expect_equal(meanTurningAngle(beaucheminTrack(1000,p.persist=1)), 0)
} )
