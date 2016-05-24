# testthat enchancement 14 (Add check on values of ufr, alpha) has been implemented
#

context("Enhancement 14 - Add check on values of ufr, alpha")

test_that("Check on value of ufr works", {
	
	TimesVector <- c(1, 2, 3, 4)
	CashflowMatrix <- diag( 4 )
	MarketValueVector <- exp( -0.04 * TimesVector )
	ufr <- -1
	alpha <- 0.1
	
	expect_that( fFitSmithWilsonYieldCurve( TimesVector, CashflowMatrix, MarketValueVector, ufr, alpha ),
				 gives_warning( "Parameter ufr should be greater than zero" ) )
				 	
})


# test_that("Check on value of ufr does not generate false positive", {
# 	
# 	TimesVector <- c(1, 2, 3, 4)
# 	CashflowMatrix <- diag( 4 )
# 	MarketValueVector <- exp( -0.04 * TimesVector )
# 	ufr <- 0.04
# 	alpha <- 0.1
# 	
# 	expect_that( fFitSmithWilsonYieldCurve( TimesVector, CashflowMatrix, MarketValueVector, ufr, alpha ),
# 				 !gives_warning( "Parameter ufr (ultimate forward rate) should be greater than zero") )
# 	
# })


test_that("Check on value of alpha works, alpha < 0 ", {
	
	TimesVector <- c(1, 2, 3, 4)
	CashflowMatrix <- diag( 4 )
	MarketValueVector <- exp( -0.04 * TimesVector )
	ufr <- 0.04
	alpha <- -0.1
	
	expect_that( fFitSmithWilsonYieldCurve( TimesVector, CashflowMatrix, MarketValueVector, ufr, alpha ),
				 throws_error( "Parameter alpha must be greater than zero") )
	
})



test_that("Check on value of alpha works, alpha = 0", {
	
	TimesVector <- c(1, 2, 3, 4)
	CashflowMatrix <- diag( 4 )
	MarketValueVector <- exp( -0.04 * TimesVector )
	ufr <- 0.04
	alpha <- 0.0
	
	expect_that( fFitSmithWilsonYieldCurve( TimesVector, CashflowMatrix, MarketValueVector, ufr, alpha ),
				 throws_error( "Parameter alpha must be greater than zero") )
	
})

