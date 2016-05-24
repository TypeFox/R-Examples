# testthat we can replicate the EIOPA (CEIOPS) QIS5 examples
#
# For the source of the tests in this file, see http://eiopa.europa.eu/fileadmin/tx_dam/files/consultations/QIS/QIS5/ceiops-paper-extrapolation-risk-free-rates_en-20100802.pdf
#


context("Replicate EIOPA QIS5 Examples")

test_that("fFitSmithWilsonYieldCurveToInstruments replicates EIOPA QIS5 example 1", {
	
	InstrumentSet1 <- read.csv("InstrumentSet1.csv")
	P4 <- 0.885
	ufr <- log( 1 + 0.042 )
	alpha <- 0.1
	reportedDigits <- 3
		
	Curve <- fFitSmithWilsonYieldCurveToInstruments(InstrumentSet1, ufr, alpha )
	
	expect_that( as.numeric( round( Curve$P( 4 ), reportedDigits ) ), equals(P4) )
	
	
})

test_that("fFitSmithWilsonYieldCurveToInstruments replicates EIOPA QIS5 example 2", {
	
	InstrumentSet1 <- read.csv("InstrumentSet2.csv")
	P4 <- 0.8836
	ufr <- log( 1 + 0.042 )
	alpha <- 0.1
	reportedDigits <- 4
	
	Curve <- fFitSmithWilsonYieldCurveToInstruments(InstrumentSet1, ufr, alpha )
	
	expect_that( as.numeric( round( Curve$P( 4 ), reportedDigits ) ), equals(P4) )
	
	
})