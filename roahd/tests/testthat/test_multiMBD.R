
# TESTING MULTIVARIATE DEPTHS ---------------------------------------------

N = 1e2
P = 1e3

time_grid = seq( 0, 10, length.out = P )

Cov = exp_cov_function( time_grid, alpha = 0.2, beta = 0.8 )

Data_1 = generate_gauss_fdata( N, centerline = rep( 0, P ), Cov = Cov )

Data_2 = generate_gauss_fdata( N, centerline = rep( 0, P ), Cov = Cov )

# Multivariate modified band depths

test_that( 'Multivariate Modified Band Depths - with management of ties',
           expect_equal( multiMBD( list( Data_1, Data_2 ), weights = 'uniform',
                                   manage_ties = TRUE ) -
                           ( MBD( Data_1, manage_ties = TRUE ) +
                               MBD( Data_2, manage_ties = TRUE ) ) / 2,
                         rep( 0, N ) )
           )

test_that( 'Multivariate Modified Band Depths - without management of ties',
           expect_equal( multiMBD( list( Data_1, Data_2 ), weights = 'uniform',
                                   manage_ties = FALSE ) -
                           ( MBD( Data_1, manage_ties = FALSE ) +
                               MBD( Data_2, manage_ties = FALSE ) ) / 2,
                         rep( 0, N ) )
)

test_that( 'Multivariate Modified Band Depths - nonuniform weights',
           expect_equal( multiMBD( list( Data_1, Data_2 ),
                                   weights = c( 1/3, 2/3 ),
                                   manage_ties = FALSE ) -
                           ( 1 / 3  * MBD( Data_1, manage_ties = FALSE ) +
                               2 / 3 * MBD( Data_2, manage_ties = FALSE ) ),
                         rep( 0, N ) )
)

test_that( 'Multivariate Modified Band Depths - wrong weights',
           expect_error( multiMBD( list( Data_1, Data_2 ),
                                   weights = c( 1/2, 1 ) ) )
)

test_that( 'Multivariate Modified Band Depths - wrong weights (bis)',
           expect_error( multiMBD( list( Data_1, Data_2 ),
                                   weights = 'unif' ) )
)


# Multivariate band depths

test_that( 'Multivariate Band Depths - with management of ties',
           expect_equal( multiBD( list( Data_1, Data_2 ),
                                   weights = 'uniform' ) -
                           ( BD( Data_1 ) +
                               BD( Data_2 ) ) / 2,
                         rep( 0, N ) )
)

test_that( 'Multivariate Band Depths - without management of ties',
           expect_equal( multiBD( list( Data_1, Data_2 ),
                                   weights = 'uniform' ) -
                           ( BD( Data_1 ) +
                               BD( Data_2 ) ) / 2,
                         rep( 0, N ) )
)

test_that( 'Multivariate Band Depths - nonuniform weights',
           expect_equal( multiBD( list( Data_1, Data_2 ),
                                   weights = c( 1/3, 2/3 ) ) -
                           ( 1 / 3  * BD( Data_1 ) +
                               2 / 3 * BD( Data_2 ) ),
                         rep( 0, N ) )
)

test_that( 'Multivariate Band Depths - wrong weights',
           expect_error( multiBD( list( Data_1, Data_2 ),
                                   weights = c( 1/2, 1 ) ) )
)

test_that( 'Multivariate Band Depths - wrong weights (bis)',
           expect_error( multiBD( list( Data_1, Data_2 ),
                                   weights = 'unif' ) )
)

