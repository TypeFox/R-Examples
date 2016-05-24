
# TESTING MBDs ------------------------------------------------------------

time_grid = seq( 0, 1, length.out = 1e2 )


Data_ref = matrix( c( 0  + sin( 2 * pi * time_grid ),
                      1  + sin( 2 * pi * time_grid ),
                      -1 + sin( 2 * pi * time_grid ) ),
                   nrow = 3, ncol = length( time_grid ), byrow = TRUE )

Data_test_1 = matrix( c( 0.6 + sin( 2 * pi * time_grid ) ),
                      nrow = 1, ncol = length( time_grid ), byrow = TRUE )

Data_test_2 = matrix( c( 0.6 + sin( 2 * pi * time_grid ) ),
                      nrow = length( time_grid ), ncol = 1, byrow = TRUE )

Data_test_3 = 0.6 + sin( 2 * pi * time_grid )

Data_test_4 = array( 0.6 + sin( 2 * pi * time_grid ), dim = length( time_grid ) )

Data_test_5 = array( 0.6 + sin( 2 * pi * time_grid ), dim = c( 1, length( time_grid ) ) )

Data_test_6 = array( 0.6 + sin( 2 * pi * time_grid ), dim = c( length( time_grid ), 1 ) )

Data_test_7 = matrix( c( 0.5  + sin( 2 * pi * time_grid ),
                         -0.5 + sin( 2 * pi * time_grid ),
                         1.1 + sin( 2 * pi * time_grid ) ),
                      nrow = 3, ncol = length( time_grid ), byrow = TRUE )

test_that( "Correctness of relative MBD (single test function in row matrix form)",
           expect_equal( MBD_relative( Data_test_1, Data_ref ), 2/ 3  ) )

test_that( "Correctness of relative MBD (single test function in column matrix form)",
           expect_equal( MBD_relative( Data_test_2, Data_ref ), 2/ 3  ) )

test_that( "Correctness of relative MBD (single test function in vector form)",
           expect_equal( MBD_relative( Data_test_3, Data_ref ), 2/ 3  ) )

test_that( "Correctness of relative MBD (single test function in 1D array form)",
           expect_equal( MBD_relative( Data_test_4, Data_ref ), 2/ 3  ) )

test_that( "Correctness of relative MBD (single test function in row-like 2D array form)",
           expect_equal( MBD_relative( Data_test_5, Data_ref ), 2/ 3  ) )

test_that( "Correctness of relative MBD (single test function in column-like 2D array form)",
           expect_equal( MBD_relative( Data_test_6, Data_ref ), 2/ 3  ) )

test_that( "Correctness of relative MBD (multiple test function)",
           expect_equal( MBD_relative( Data_test_7, Data_ref ), c( 2/3, 2/3, 0 ) ))




