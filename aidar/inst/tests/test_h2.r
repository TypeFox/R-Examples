
histoFile = system.file("extdata", "histos.xml.gz", package="aidar")

test_that("reading a 2D histogram from file", {

	h10 = getHisto2D(histoFile, '10')
	
	expect_equal( class(h10), "data.frame" )
	expect_equal( typeof(h10), "list" )
	expect_equal( length(h10), 9 )
	expect_equal( names(h10) , c("binNumberX", "binNumberY", "binX", "binY", "entries", "error", "height", "weightedMeanX", "weightedMeanY") )
	
	expect_equal( mean(h10$height)      , mean(h10$entries) )
	
	expect_equal( mean(h10$entries)      ,  8.77193, tolerance = 1.E-6, scale = mean(h10$entries) )
	expect_equal( mean(h10$weightedMeanX), 25.40455,  tolerance = 1.E-6, scale = mean(h10$weightedMeanX) )
	expect_equal( mean(h10$weightedMeanY), 50.14967, tolerance = 1.E-6, scale = mean(h10$weightedMeanY) )
	expect_equal( mean(h10$error)        ,  2.53375, tolerance = 1.E-6, scale = mean(h10$error)  )

})


test_that("getAnnotation(2D histo)", {

	a10 = getAnnotation(histoFile, '10')

	expect_equal( length(a10$key), 8 )
	
	expect_equal( a10$values[[3]], "993" )
	expect_equal( a10$values[[4]], "24.985" )	
})
