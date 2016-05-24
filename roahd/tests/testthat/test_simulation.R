


# TESTING THE SIMULATION OF UNIVARIATE FUNCTIONAL DATA --------------------

N = 30
P = 1e2

t0 = 0
tP = 1

time_grid = seq( t0, tP, length.out = P )

C = exp_cov_function( time_grid, alpha = 0.1, beta = 0.2 )

CholC = chol( C )

centerline = sin( 2 * pi * time_grid )

test_that( 'Generation of gaussian univariate functional data - 1',
           expect_silent( generate_gauss_fdata( N,
                                                 centerline,
                                                 Cov = C )
           ) )


test_that( 'Generation of gaussian univariate functional data - 1',
           expect_silent( generate_gauss_fdata( N,
                                                centerline,
                                                CholCov = CholC )
           ) )



# TESTING THE SIMULATION OF MULTIVARIATE FUNCTIONAL DATA ------------------

N = 30
P = 1e2
L = 3

t0 = 0
tP = 1

time_grid = seq( t0, tP, length.out = P )

C1 = exp_cov_function( time_grid, alpha = 0.1, beta = 0.2 )
C2 = exp_cov_function( time_grid, alpha = 0.2, beta = 0.5 )
C3 = exp_cov_function( time_grid, alpha = 0.3, beta = 1 )

CholC1 = chol( C1 )
CholC2 = chol( C2 )
CholC3 = chol( C3 )

centerline = matrix( c( sin( 2 * pi * time_grid ),
                        sqrt( time_grid ),
                        10 * ( time_grid - 0.5 ) * time_grid ),
                     nrow = 3, byrow = TRUE )

test_that( 'Generation of gaussian multivarite functional data - 1',
           expect_silent( generate_gauss_mfdata( N, L,
                                                 centerline,
                                                 correlations = c( 0.5, 0.5,
                                                                   0.5 ),
                                                 listCov = list( C1, C2, C3 ) )
                          ) )

test_that( 'Generation of gaussian multivarite functional data - 2',
           expect_silent( generate_gauss_mfdata( N, L,
                                                 centerline,
                                                 correlations = c( 0.5, 0.5,
                                                                   0.5 ),
                                                 listCholCov = list( CholC1,
                                                                     CholC2,
                                                                     CholC3 ) )
                          ) )

test_that( 'Generation of gaussian multivarite functional data - 3',
           expect_error( generate_gauss_mfdata( N, L,
                                                 centerline,
                                                 correlations = c( 0.5, 0.5,
                                                                   0.5 ),
                                                 listCholCov = list( CholC1[-1,],
                                                                     CholC2[-1,],
                                                                     CholC3[-1,]
                                                                     ) ) ) )

test_that( 'Generation of gaussian multivarite functional data - 4',
           expect_error( generate_gauss_mfdata( N, L,
                                                centerline[-1,],
                                                correlations = c( 0.5, 0.5,
                                                                  0.5 ),
                                                listCov = c( C1, C2, C3 ),
                                                listCholCov = list( CholC1[-1,],
                                                                    CholC2[-1,],
                                                                    CholC3[-1,]
                                                                    ) ) ) )

test_that( 'Generation of gaussian multivarite functional data - 5',
           expect_error( generate_gauss_mfdata( N, L,
                                                centerline[,-1],
                                                correlations = c( 0.5, 0.5 ),
                                                listCov = c( C1, C2, C3 ),
                                                listCholCov = list( CholC1,
                                                                    CholC2,
                                                                    CholC3 ) )
                         ) )
