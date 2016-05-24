
histoFile = system.file("extdata", "clouds.xml.gz", package="aidar")

test_that("reading a 1D cloud from file", {

	c21 = getCloud1D(histoFile, '21')
	
	expect_equal( class(c21), "data.frame" )
	expect_equal( typeof(c21), "list" )
	expect_equal( length(c21), 1 )
	expect_equal( names(c21) , c("valuesX") )
	expect_equal( length(c21$valuesX), 100 )
	
	expect_equal( mean(c21$valuesX), 26.12283, tolerance = 1.E-6, scale = mean(c21$valuesX) )
	
})


test_that("getAnnotation(1D cloud)", {

	ann = getAnnotation(histoFile, '21')

	expect_equal( length(ann$keys), 5 )
	
	expect_equal( ann$values[ann$keys=="Entries"], "100" )
	expect_equal( ann$values[ann$keys=="Mean"]   , "26.123" )
})
