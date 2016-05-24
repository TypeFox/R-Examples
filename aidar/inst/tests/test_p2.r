
histoFile = system.file("extdata", "histos.xml.gz", package="aidar")

test_that("reading a 2D profile histogram from file", {

	p2d = getProfile2D(histoFile, 'Example 2D profile (gauss)')
	
	expect_equal( class(p2d), "data.frame" )
	expect_equal( typeof(p2d), "list" )
	expect_equal( length(p2d), 10 )
	expect_equal( names(p2d) , c("binNumberX", "binNumberY", "binX", "binY", "entries", "error", "height", "rms", "weightedMeanX", "weightedMeanY") )
	
	expect_equal( mean(p2d$entries)      , 9.259259, tolerance = 1.E-6, scale = mean(p2d$entries) )
	expect_equal( mean(p2d$height)       , 48.76922, tolerance = 1.E-6, scale = mean(p2d$height) )
	expect_equal( mean(p2d$weightedMeanX), 25.27984, tolerance = 1.E-6, scale = mean(p2d$weightedMeanX) )
	expect_equal( mean(p2d$weightedMeanY), 21.33456, tolerance = 1.E-6, scale = mean(p2d$weightedMeanY) )
	expect_equal( mean(p2d$rms)          , 7.943256, tolerance = 1.E-6, scale = mean(p2d$rms)  )
	expect_equal( mean(p2d$error)        , 3.748361, tolerance = 1.E-6, scale = mean(p2d$error)  )

})

test_that("getAnnotation(2D profile)", {

	ann = getAnnotation(histoFile, 'Example 2D profile (gauss)')

	expect_equal( length(ann$key), 8 )
	
	expect_equal( ann$values[[3]], "953" )
	expect_equal( ann$values[[4]], "24.748" )	
})
