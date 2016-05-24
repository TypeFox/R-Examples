
histoFile = system.file("extdata", "clouds.xml.gz", package="aidar")

test_that("reading a 3D cloud from file", {

	c3d = getCloud3D(histoFile, '33')
	
	expect_equal( class(c3d), "data.frame" )
	expect_equal( typeof(c3d), "list" )
	expect_equal( length(c3d), 3 )
	expect_equal( names(c3d) , c("valuesX", "valuesY", "valuesZ") )
	expect_equal( length(c3d$valuesX), 100 )
	expect_equal( length(c3d$valuesX), length(c3d$valuesY) )
	expect_equal( length(c3d$valuesX), length(c3d$valuesZ) )
	
	expect_equal( mean(c3d$valuesX), 26.12283, tolerance = 1.E-6, scale = mean(c3d$valuesX) )
	expect_equal( mean(c3d$valuesY), 50.28191, tolerance = 1.E-6, scale = mean(c3d$valuesY) )
	expect_equal( mean(c3d$valuesZ), 0.1501225, tolerance = 1.E-6, scale = mean(c3d$valuesZ) )
	
})


test_that("getAnnotation(3D cloud)", {

	ann = getAnnotation(histoFile, '33')

	expect_equal( length(ann$keys), 9 )

	expect_equal( ann$values[ann$keys=="Entries"], "100" )
	expect_equal( ann$values[ann$keys=="MeanX"]  , "26.123" )
	expect_equal( ann$values[ann$keys=="MeanY"]  , "50.282" )
	expect_equal( ann$values[ann$keys=="MeanZ"]  , "0.1501" )

})

