
histoFile = system.file("extdata", "histos.xml.gz", package="aidar")

test_that("getFileInfo", {

	info = getFileInfo(histoFile)

	expect_equal( length(info$histogram1d$name), 2 )

	expect_equal( info$histogram1d$name   [[2]], "2" )
	expect_equal( info$histogram1d$title  [[2]], "Example histogram 2 (flat)" )

	expect_equal( info$histogram1d$entries[[1]], "984" )
	
})

