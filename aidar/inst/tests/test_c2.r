
histoFile = system.file("extdata", "clouds.xml.gz", package="aidar")

test_that("reading a 2D cloud from file", {

	c2d = getCloud2D(histoFile, '30')
	
	expect_equal( class(c2d), "data.frame" )
	expect_equal( typeof(c2d), "list" )
	expect_equal( length(c2d), 2 )
	expect_equal( names(c2d) , c("valuesX", "valuesY") )
	expect_equal( length(c2d$valuesX), 100 )
	expect_equal( length(c2d$valuesX), length(c2d$valuesY) )
	
	expect_equal( mean(c2d$valuesX), 26.12283, tolerance = 1.E-6, scale = mean(c2d$valuesX) )
	expect_equal( mean(c2d$valuesY), 50.28191, tolerance = 1.E-6, scale = mean(c2d$valuesY) )
	
})


test_that("getAnnotation(2D cloud)", {

	ann = getAnnotation(histoFile, '30')

	expect_equal( length(ann$keys), 7 )
	
	expect_equal( ann$values[ann$keys=="Entries"], "100" )
	expect_equal( ann$values[ann$keys=="MeanX"]  , "26.123" )
	expect_equal( ann$values[ann$keys=="MeanY"]  , "50.282" )
})
