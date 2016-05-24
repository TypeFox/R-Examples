
histoFile = system.file("extdata", "histos.xml.gz", package="aidar")

test_that("reading a 1D histogram from file", {

	h1 = getHisto1D(histoFile, '1')
	
	expect_equal( class(h1), "data.frame" )
	expect_equal( typeof(h1), "list" )
	expect_equal( length(h1), 6 )
	expect_equal( names(h1) , c("binNumber","binX","entries","error","height","weightedMean") )
	
	expect_equal( mean(h1$height)      , mean(h1$entries) )
	
	expect_equal( mean(h1$entries)     , 19.60784, tolerance = 1.E-6, scale = mean(h1$entries) )
	expect_equal( mean(h1$weightedMean), 25.54948, tolerance = 1.E-6, scale = mean(h1$weightedMean) )
	expect_equal( mean(h1$error)       , 4.103032, tolerance = 1.E-6, scale = mean(h1$error)  )

})

test_that("getAnnotation(1D histo)", {

	a10 = getAnnotation(histoFile, '1')

	expect_equal( length(a10$key), 8 )
	
	expect_equal( a10$values[[3]], "984" )
	expect_equal( a10$values[[4]], "24.705" )	
})
