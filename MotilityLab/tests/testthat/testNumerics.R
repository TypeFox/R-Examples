test_that("Potential numerical issues",{
	expect_equal( range( sapply(subtracks( diffinv(matrix(rep(sqrt(2),200),ncol=4,byrow=T)), 2 ),
		overallAngle ) ), c(0,0) )
})