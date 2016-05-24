
tupleFile = system.file("extdata", "tuple.xml.gz", package="aidar")

test_that("reading a tuple from file", {

	t100 = getTuple(tupleFile, '100')
	
	expect_equal( class(t100), "data.frame" )
	expect_equal( typeof(t100), "list" )
	expect_equal( names(t100) , c( "lin", "sin", "gauss", "flat" ) )

	expect_equal( length(t100),    4 )  # nCol
	expect_equal( length(t100$lin), 1000 )  # nRow	
	
	expect_equal( mean(t100$sin)  ,  0.001628846, tolerance = 1.E-6, scale = mean(t100$sin) )
	expect_equal( mean(t100$gauss), -0.000790439, tolerance = 1.E-6, scale = mean(t100$gauss) )
	expect_equal( mean(t100$flat) , 24.89807    , tolerance = 1.E-6, scale = mean(t100$flat) )
	expect_equal( max (t100$lin)  , 999 )
	
})

test_that("getAnnotation", {

	a100 = getAnnotation(tupleFile, '100')

	expect_equal( length(a100), 2 )

	expect_equal( a100$keys[[1]]  , "Title" )
	expect_equal( a100$values[[1]], "100" )
	
})

