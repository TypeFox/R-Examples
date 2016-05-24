test.computeNextWeight = function()
{
	checkEquals( computeNextWeight( 1, 1 ), c( 0.5, 0.5 ), " equal")
	checkEquals( computeNextWeight( 0, 1 ), c( 1, 0 ), " one 0")
	checkEquals( computeNextWeight( 4, 4/3 ), c( 0.25, 0.75 ), " mixed")

}

test.computeModel = function()
{
	# unfortunately everything done in 1 function is now necessary for speed with fitting, but very hard to test
	ones = TimedVector( c( 1, 1, 1, 1, 1 ), c( 1, 2, 3, 63, 64 ) )

	moe = computeModel( ones, mFast = 1, mSlow = 1, n = 1, g = 0,
		h = 0.5, tau = 1/60, threshold = 10 )
	est = attr( moe, "estimates" )
	mse = attr( moe, "errorEstimates" )
	weights = attr( moe, "weights" )
	y = attr( moe, "y" )
	
	checkEquals( est, cbind( c( 0, 1, 1, 1, 1 ), c( 0, 1, 1, 1, 1 ) ), " estimates" )
	checkEquals( mse, cbind( c( 1, 1, 0, 0, 0 ), c( 1, 1, 0, 0, 0 ) ), " errors" )
	checkEquals( weights, cbind( rep( 0.5, 5 ), rep( 0.5, 5 ) ), " weights" )
	checkEquals( y, rep( 1, 5 ), " y" )
	checkEquals( as.vector( moe ), c( 0, 1, 1, 1, 1 ), " forecast" )
	checkEquals( time( moe ), c( 1, 2, 3, 63, 64 ), " time" )

	moe = computeModel( ones, mFast = 0.8, mSlow = 0.02, n = 0.05, g = 10,
		h = 0.1, tau = 1/60, threshold = 10 )
	est = attr( moe, "estimates" )
	mse = attr( moe, "errorEstimates" )
	weights = attr( moe, "weights" )
	y = attr( moe, "y" )
	
	checkEquals( est, cbind( c( 0, 0.0000800000000000, 0.0728014545454546, 0.1964281384644158, 0.3449199231925152 ), 
							 c( 0, 0.00000200000000000, 0.00182017818181818, 0.00514744993852983, 0.00974341197289068 ) ), " estimates" )
	checkEquals( mse, cbind( c( 1, 1, 0.999992000320000, 0.992977257438651, 0.975614781399307 ), 
							 c( 1, 1, 0.999999800000200, 0.099981795783444, 0.144469285812462 ) ), " errors" )
	checkEquals( weights, cbind( c( 0.5, 0.5, 0.500001949928044, 0.091478080069599, 0.128980752464492 ), 
								  c( 0.5, 0.5, 0.499998050071956, 0.908521919930401, 0.871019247535508 ) ), " weights" )
	checkEquals( y, c( 0.000100000000000, 0.090909090909091, 0.166666951384122, 0.230987096232345, 0.286195436224904 ), " y" )
	checkEquals( as.vector( moe ), c( 0, 0.0000410000000000, 0.0373109547720178, 0.0226454400792688, 0.0529747305984210 ), " forecast" )
	checkEquals( time( moe ), c( 1, 2, 3, 63, 64 ), " time" )
	
}

test.averageBySession = function()
{
	responses = c( 3, 2, 5, 2, 3, 2, 3, 1, 1 )
	sessions = c( 1, 6, 8, 10 )
	
	checkEquals( averageBySession( responses, sessions ), c( 3, 2.5, 1 ), " average responses" )
}

test.averageBySessionNA = function()
{
	responses = c( 3, 2, 5, 2, 3, 2, 3, 1, 1 )
	sessions = NA
	
	checkEquals( averageBySession( responses, sessions ), responses, " same responses" )
}

test.calculateResponse = function()
{
	est = c(0, 0.5, 1)
	checkEquals( calculateResponse( 0, 10, est ), c( NaN, 10, 10 ), " k = 0" )
	checkEquals( calculateResponse( 0.1, 10, est ), c( 0, 9, 10 ), " k = 0.1", tolerance = 0.01 )
	checkEquals( calculateResponse( 0.5, 10, est ), c( 0, 10 * 2 / 3, 10 ), " k = 0.5" )
	checkEquals( calculateResponse( 1, 10, est ), est * 10, " k = 1" )
}

test.fitModel = function()
{
	trials = TimedVector(c(rep(1, 20), rep(0, 20)), c(1:20, (2*TV_DAY):(2*TV_DAY + 4), (4*TV_DAY):(4*TV_DAY + 4),
													   (6*TV_DAY):(6*TV_DAY + 4), (8*TV_DAY):(8*TV_DAY + 4)))
	params = c(0.9, 0.2, 0.1, 0.1, 1, 10, 0)
	est = computeModel(trials, mFast=params[1], mSlow=params[2], n=params[3], h=params[4], g=params[7], verbose = TRUE)
	responses = calculateResponse(params[5], params[6], est)
	# don't estimate g since it's a different scale thatn the other params
	results = fitModel(dataX = list(trials), dataResponse = list(responses), fitG = FALSE )
	fitParams = results$par[1,]
	# large tolerance since we're just checking running fit and can vary due to stochastic nature of algorithm
	checkEqualsNumeric( fitParams, params[1:6], " params compare", tolerance = 0.5 )
}

computeNextWeight = function( fastMse, slowMse )
{
	value = .C( "compute_weight",
				as.double( fastMse ),
				as.double( slowMse ),
				fastWeight = double( 1 ),
				slowWeight = double( 1 ) )
	return( c( value$fastWeight, value$slowWeight ) )
}
