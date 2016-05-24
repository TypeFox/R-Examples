
# TESTING OUTLIERGRAM -----------------------------------------------------

set.seed( 1618 )

N = 200
N_outliers = 4
P = 200

grid = seq( 0, 1, length.out = P )


Cov = exp_cov_function( grid, alpha = 0.2, beta = 0.8 )

Data = generate_gauss_fdata( N,
                             centerline = sin( 4 * pi * grid ),
                             Cov = Cov )


Data_out = array( 0, dim = c( N_outliers, P ) )

Data_out[ 1, ] = generate_gauss_fdata( 1,
                                       sin( 4 * pi * grid + pi / 2 ),
                                       Cov = Cov )

Data_out[ 2, ] = generate_gauss_fdata( 1,
                                       sin( 4 * pi * grid - pi / 2 ),
                                       Cov = Cov )

Data_out[ 3, ] = generate_gauss_fdata( 1,
                                       sin( 4 * pi * grid + pi/ 3 ),
                                       Cov = Cov )

Data_out[ 4, ] = generate_gauss_fdata( 1,
                                       sin( 4 * pi * grid - pi / 3),
                                       Cov = Cov )
Data = rbind( Data, Data_out )

fD = fData( grid, Data )


test_that( 'Outliergram - general',
           expect_silent( outliergram( fD, display = FALSE ) ) )

test_that( 'Outliergram - no adjustment',
           expect_equal(
             outliergram( fD, display = FALSE )$ID_outliers
             ,
             c( 31,  78, 117, 122, 152, 183, 201, 202, 203, 204 ) ) )

test_that( 'Outliergram - no adjustment 2',
           expect_equal(
             outliergram( fD, Fvalue = 10,
                          display = FALSE )$ID_outliers,
             c( 201, 202, 203, 204 ) ) )


test_that( 'Outliergram - with adjustment',{
  testthat::skip_on_cran()
           expect_equal(
             outliergram( fD,
                          adjust = list( N_trials = 10,
                                         trial_size = 5 * nrow( Data ),
                                         TPR = 0.01,
                                         VERBOSE = FALSE ),
                          display = FALSE )$ID_outliers,
             c( 78, 117, 183, 201, 202, 203, 204 ) ) })

