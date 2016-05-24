# Test that defect 17 (Error in creating cashflow matrix) has been fixed
# Note this is a user / input error, we now include a meaningful error message

context( "Defect 17 (Error in creating cashflow matrix)" )

test_that( "Defect 17 is fixed", {
	
	expect_that( fGetCashflowsSwap( data.frame( Frequency=1, Rate=4, Tenor=0.5 ) ),
				 throws_error( "No cashflows calculated for swap, check Tenor and Frequency" ) ) 
		
})