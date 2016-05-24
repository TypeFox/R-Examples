library('testthat')
context("createProductRuleGrid")

test_that("check for sensible input values", {
    correct_type1 <- "GQU"
    correct_dimension <- 3
    correct_k <- 4
    correct_sym <- TRUE
    
    wrong_type1 <- "HGE"
    wrong_type2 <- function() { return( list( "nodes"=c(1), "weights"=c(1) ) ) }			# does not contain argument
    wrong_type3 <- function(n1, n2) { return( list( "nodes"=c(1), "weights"=c(1) ) ) }	    # contains too many arguments
    wrong_type4 <- function(n) { return( c(1) ) }											# does not return list
    wrong_type5 <- function(n) { return( list( "weights"=c(1) ) ) }							# return list does not contain "nodes"
    wrong_type6 <- function(n) { return( list( "nodes"=c(1) ) ) }						    # return list does not contain "weights"
    
	wrong_dimension <- 3.1
	wrong_k         <- 4.1
	wrong_sym1      <- NA
	wrong_sym2      <- 1
	
    expect_that(createProductRuleGrid( type=wrong_type1,   dimension=correct_dimension, k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createProductRuleGrid( type=wrong_type2,   dimension=correct_dimension, k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createProductRuleGrid( type=wrong_type3,   dimension=correct_dimension, k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createProductRuleGrid( type=wrong_type4,   dimension=correct_dimension, k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createProductRuleGrid( type=wrong_type5,   dimension=correct_dimension, k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createProductRuleGrid( type=wrong_type6,   dimension=correct_dimension, k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createProductRuleGrid( type=correct_type1, dimension=wrong_dimension,   k=correct_k, sym=correct_sym ), throws_error())
    expect_that(createProductRuleGrid( type=correct_type1, dimension=correct_dimension, k=wrong_k,   sym=correct_sym ), throws_error())
    expect_that(createProductRuleGrid( type=correct_type1, dimension=correct_dimension, k=correct_k, sym=wrong_sym1 ),  throws_error())
    expect_that(createProductRuleGrid( type=correct_type1, dimension=correct_dimension, k=correct_k, sym=wrong_sym2 ),  throws_error())
  
})

test_that("check dimensions of grid", {
	
	type            <- "GQU"
	accuracy		<- 5

	intgrid 		<- createProductRuleGrid( type, dimension=1, k=5 )
	num.1d.nodes	<- length( intgrid$weights )
	for (dimension in 1:4) {
		intgrid <- createProductRuleGrid( type, dimension=dimension, k=accuracy )
		expect_that( nrow( intgrid$nodes ),     equals( num.1d.nodes^dimension ) )
		expect_that( ncol( intgrid$nodes ),     equals( dimension ) )
		expect_that( length( intgrid$weights ), equals( num.1d.nodes^dimension ) )
	}
})
