library('testthat')
context("sortrows")

test_that("check for sensible input values", {
    A_vector <- c(1,1,1)
    A_list   <- list( c(1,2,3), c(2,3,4) )
    
    expect_that(sortrows( A = A_vector ), throws_error())
    expect_that(sortrows( A = A_list ), throws_error())
  
})

test_that("check output values for some matrices", {

    # test with 1x3 matrix
    A1          <- matrix( c(1,3,2), nrow=1, ncol=3, byrow=TRUE )
    A1.sorted   <- A1
    expect_that( sortrows( A1 ), equals( A1.sorted ) )
    
    # test with 3x1 matrix
    A2          <- matrix( c(1,3,2), nrow=3, ncol=1, byrow=TRUE )
    A2.sorted   <- matrix( c(1,2,3), nrow=3, ncol=1, byrow=TRUE )
    expect_that( sortrows( A2 ), equals( A2.sorted ) )
    
    # test with 6x2 matrix
    A3 <- matrix( c( 
            1, 2,
            2, 2,
            3, 2,
            3, 1,
            1, 1,
            2, 1 ), nrow=6, ncol=2, byrow=TRUE )
    
    A3.sorted <- matrix( c( 
            1, 1,
            1, 2,
            2, 1,
            2, 2,
            3, 1,
            3, 2 ), nrow=6, ncol=2, byrow=TRUE )	
	
    expect_that( sortrows( A3 ), equals( A3.sorted ) )
	
})

test_that("check returned indices", {

    # test with 6x2 matrix
    A3 <- matrix( c( 
            1, 2,
            2, 2,
            3, 2,
            3, 1,
            1, 1,
            2, 1 ), nrow=6, ncol=2, byrow=TRUE )
    
    A3.sorted <- matrix( c( 
            1, 1,
            1, 2,
            2, 1,
            2, 2,
            3, 1,
            3, 2 ), nrow=6, ncol=2, byrow=TRUE )	
	
    # calculate result including vector of indices
    res <- sortrows( A3, index.return=TRUE )
    
    expect_that( res, is_a( "list" ) )                  # res is a list
    expect_that( "x" %in% names(res), is_true() )       # res$x exists
    expect_that( "ix" %in% names(res), is_true() )      # res$ix exists

    expect_that( A3[ res$ix, ], equals( A3.sorted ) )    # indices can be used to sort rows
})

