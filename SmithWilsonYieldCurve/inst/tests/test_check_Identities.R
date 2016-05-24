# Test that when using one calibration instrument that we retrieve the same price

context("Check Identities")

test_that("Single calibration rate is retrieved", {

	testRate <- 0.05
	testZCB <- exp( -log( 1 + testRate ) )
	
	InstrumentSet <- data.frame( Type="LIBOR", Tenor=1, Rate=testRate, Frequency=NA )
	
	Curve <- fFitSmithWilsonYieldCurveToInstruments(InstrumentSet, 0.04, 0.1)
	
	expect_that( as.numeric( round( Curve$P(1), digits=5 ) ), equals( round( testZCB, digits=5 ) ) )
	
} )