
histoFile = system.file("extdata", "histos.xml.gz", package="aidar")

test_that("reading a 3D histogram from file", {

	h3d = getHisto3D(histoFile, '13')
	
	expect_equal( class(h3d), "data.frame" )
	expect_equal( typeof(h3d), "list" )
	expect_equal( length(h3d), 12 )
	expect_equal( names(h3d) , c("binNumberX", "binNumberY", "binNumberZ", "binX", "binY", "binZ", "entries", "error", "height", "weightedMeanX", "weightedMeanY", "weightedMeanZ") )
	
	expect_equal( mean(h3d$height)      , mean(h3d$entries) )
	
	expect_equal( mean(h3d$entries)      , 7.246377, tolerance = 1.E-6, scale = mean(h3d$entries) )
	expect_equal( mean(h3d$weightedMeanX), 25.14937,  tolerance = 1.E-6, scale = mean(h3d$weightedMeanX) )
	expect_equal( mean(h3d$weightedMeanY), 50.51821, tolerance = 1.E-6, scale = mean(h3d$weightedMeanY) )
	expect_equal( mean(h3d$weightedMeanZ), 29.48783, tolerance = 1.E-6, scale = mean(h3d$weightedMeanZ) )
	expect_equal( mean(h3d$error)        , 2.269346, tolerance = 1.E-6, scale = mean(h3d$error)  )

})

test_that("getAnnotation(3D histo)", {

	ann = getAnnotation(histoFile, '13')

	expect_equal( length(ann$key), 10 )
	
	expect_equal( ann$values[[3]], "992" )
	expect_equal( ann$values[[4]], "24.986" )	

})

