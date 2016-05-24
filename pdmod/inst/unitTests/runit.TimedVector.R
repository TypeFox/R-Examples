test.create = function()
{
	checkEquals( TimedVector( c( 1, 1, 1) ), TimedVector( c( 1, 1, 1 ), 1:3 ), " missing time")
	checkEquals( as.vector( TimedVector( c( 3, 1, 5), 1:3 ) ), c( 3, 1, 5), " values")
	checkEquals( TimedVector( NULL, NULL ), NULL, " null")
	
	#errors
	checkException( TimedVector( 1, 1:2 ), " different lengths")
	checkException( TimedVector( 1:3, 3:1 ), "decreasing time" )
	checkException( TimedVector( rep( 1, 3 ), rep( 2, 3 ) ), " duplicate time" ) # should this be an error?
}

test.time = function()
{
	checkEquals( time( TimedVector( c( 3, 1, 5), 1:3 ) ), 1:3, " time")
}

test.is = function()
{
	checkEquals( isTimedVector( TimedVector( 1 ) ), TRUE, " object" )
	checkEquals( isTimedVector( NULL ), FALSE, " null" )
}

test.c = function()
{
	a = TimedVector( 1:2, 1:2 )
	b = TimedVector( c( 1, 1, 1 ), 5:7 )
	
	checkEquals( c( a, b ), TimedVector( c( 1, 2, 1, 1, 1 ), c( 1, 2, 5, 6, 7 ) ), " test" )
	checkException( c( b, a ), " combined result not valid")
	
}

test.subset = function()
{
	# the expected shoul be a timed vector if the subset function is re-enabled
	a = TimedVector( 1:5, 11:15 )
	checkEquals( a[1:5], 1:5, " entire vector")
	checkEquals( a[1], c( 1 ), " first element")
	checkEquals( a[5], c( 5 ), " last element")
	checkEquals( a[2:4], 2:4, " middle")
}
