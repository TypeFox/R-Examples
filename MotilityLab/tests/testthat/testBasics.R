
test_that("turning angles work",{
	expect_equal( overallAngle(
		rbind(c(0,0,0,0),c(0,1,1,1))), 0 )
	expect_equal( aggregate(TCells, overallAngle, subtrack.length=1)[1,2], 0 )
	expect_equal( aggregate(TCells, overallDot, subtrack.length=1)[1,2],
		2*mean(rowSums(do.call(rbind,normalizeTracks(subtracks(TCells,1)))[,-1]^2)) )
})

test_that("maxTrackLength works",{
	expect_equal( maxTrackLength(TCells), 39 )
	expect_equal( maxTrackLength(BCells), 39 )
	expect_equal( maxTrackLength(Neutrophils), 55 )
})