
histoFile = system.file("extdata", "histos.xml.gz", package="aidar")

test_that("reading a 1D profile histogram from file", {

	h10 = getProfile1D(histoFile, 'Example profile (gauss)')
	
	expect_equal( class(h10), "data.frame" )
	expect_equal( typeof(h10), "list" )
	expect_equal( length(h10), 7 )
	expect_equal( names(h10) , c("binNumber", "binX", "entries", "error", "height", "rms", "weightedMean") )
	
	expect_equal( mean(h10$entries)     , 19.60784, tolerance = 1.E-6, scale = mean(h10$entries) )
	expect_equal( mean(h10$height)      , 30.18176, tolerance = 1.E-6, scale = mean(h10$height) )
	expect_equal( mean(h10$weightedMean), 24.64177, tolerance = 1.E-6, scale = mean(h10$weightedMean) )
	expect_equal( mean(h10$rms)         , 9.439369, tolerance = 1.E-6, scale = mean(h10$rms)  )
	expect_equal( mean(h10$error)       , 2.974591, tolerance = 1.E-6, scale = mean(h10$error)  )

})

test_that("getAnnotation(1D profile)", {

	a10 = getAnnotation(histoFile, 'Example profile (gauss)')

	expect_equal( length(a10$key), 8 )
	
	expect_equal( a10$values[[3]], "971" )
	expect_equal( a10$values[[4]], "20.288" )	
})
