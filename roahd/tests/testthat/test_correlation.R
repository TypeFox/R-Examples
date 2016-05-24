
# TESTING MAX-MIN FUNCTIONS -----------------------------------------------

P = 1e4

time_grid = seq( 0, 1, length.out = P )

h = time_grid[ 2 ] - time_grid[ 1 ]

Data = matrix( c( 1 * time_grid,
                  2 *  time_grid,
                  3 * ( 0.5 - abs( time_grid - 0.5 ) ) ),
               nrow = 3, ncol = P, byrow = TRUE )

fD = fData( time_grid, Data )

test_that( 'Max function for functional data, which = TRUE, grid',
           expect_equal( maxima( fD, which = TRUE )$grid,
                             c( 1, 1, 0 + 4999 * h ) ) )

test_that( 'Max function for functional data, which = TRUE, value',
           expect_equal( maxima( fD, which = TRUE )$value,
                             c( 1, 2, 3 * ( 0.5 - abs( 0.5 - 4999 * h) ) ) ) )

test_that( 'Min function for functional data, which = TRUE, grid',
           expect_equal( minima( fD, which = TRUE )$grid,
                         c( 0, 0, 0 ) ) )

test_that( 'Min function for functional data, which = TRUE, value',
           expect_equal( minima( fD, which = TRUE )$value,
                         c( 0, 0, 0 ) ) )

test_that( 'Max function for functional data, which = FALSE',
           expect_equal( maxima( fD, which = FALSE ),
                         c( 1, 2, 3 * ( 0.5 - abs( 0.5 - 4999 * h) ) ) ) )

test_that( 'Min function for functional data, which = FALSE',
           expect_equal( minima( fD, which = FALSE ),
                         c( 0, 0, 0 ) ) )

test_that( 'Area under the curve - 1',
           expect_equal( area_under_curve( fD ),
                         c( 0.5, 1, 0.75 ) ) )

fD = fData( time_grid,
            matrix( c( sin( 2 * pi * time_grid ),
                       cos( 2 * pi * time_grid ),
                       4 * time_grid * ( 1 - time_grid ) ),
                    nrow = 3, ncol = P, byrow = TRUE ) )

test_that( 'Area under the curve - 2',
           expect_true(
             all( c( area_under_curve( fD )[1:2],
                     abs( area_under_curve( fD[ 3, ] ) - 2/3 ) ) <=
                    .Machine$double.eps^0.5 ) ) )


# ORDERING ----------------------------------------------------------------

P = 1e3

time_grid = seq( 0, 1, length.out = P )

h = time_grid[ 2 ] - time_grid[ 1 ]

Data_1 = matrix( c( 1 * time_grid,
                    2 *  time_grid ),
                 nrow = 2, ncol = P, byrow = TRUE )

Data_2 = matrix( 3 * ( 0.5 - abs( time_grid - 0.5 ) ),
                 nrow = 1, byrow = TRUE )

Data_3 = rbind( Data_1, Data_1 )


fD_1 = fData( time_grid, Data_1 )
fD_2 = fData( time_grid, Data_2 )
fD_3 = fData( time_grid, Data_3 )

# MAX ORDERING
test_that( 'Max_ordering - case 1',
           expect_equal( max_ordered( fD_1, fD_2 ),
                         c( TRUE, FALSE ) ) )

test_that( 'Max_ordering - case 2',
           expect_equal( max_ordered( fD_2, fD_1 ),
                         c( FALSE, TRUE ) ) )

test_that( 'Max_ordering - case 3',
           expect_error( max_ordered( fD_1, fD_3 ) ) )

test_that( 'Max_ordering - case 4',
           expect_error( max_ordered( fD_3, fD_1 ) ) )

test_that( 'Max_ordering - case 5',
           expect_equal( max_ordered( fD_2, fD_3 ),
                         c( FALSE, TRUE, FALSE, TRUE ) ) )

test_that( 'Max_ordering - case 6',
           expect_equal( max_ordered( fD_3, fD_2 ),
                         c( TRUE, FALSE, TRUE, FALSE ) ) )

# AREA ORDERING
test_that( 'Area ordering - case 1',
           expect_equal( area_ordered( fD_1, fD_2 ),
                         c( TRUE, FALSE ) ) )

test_that( 'Area ordering - case 2',
           expect_equal( area_ordered( fD_2, fD_1 ),
                         c( FALSE, TRUE ) ) )

test_that( 'Area ordering - case 3',
           expect_error( area_ordered( fD_1, fD_3 ) ) )

test_that( 'Area ordering - case 4',
           expect_error( area_ordered( fD_3, fD_1 ) ) )

test_that( 'Area ordering - case 5',
           expect_equal( area_ordered( fD_2, fD_3 ),
                         c( FALSE, TRUE, FALSE, TRUE ) ) )

test_that( 'Area ordering - case 6',
           expect_equal( area_ordered( fD_3, fD_2 ),
                         c( TRUE, FALSE, TRUE, FALSE ) ) )


# KENDALL CORRELATION -----------------------------------------------------

N = 2e2

P = 1e3

t0 = 0
t1 = 1

time_grid = seq( t0, t1, length.out = P )

Cov = exp_cov_function( time_grid, alpha = 0.3, beta = 0.4 )

Data_1 = generate_gauss_fdata( N, centerline = sin( 2 * pi * time_grid ), Cov = Cov )
Data_2 = generate_gauss_fdata( N, centerline = sin( 2 * pi * time_grid ), Cov = Cov )

mfD = mfData( time_grid, list( Data_1, Data_2 ) )

test_that( 'Kendall correlation with max ordering ',
           expect_silent( invisible( cor_kendall( mfD, ordering = 'max' ) ) ) )

test_that( 'Kendall correlation with area ordering',
           expect_silent( invisible( cor_kendall( mfD, ordering = 'area' ) ) ) )

# SPEARMAN RANK CORRELATION -----------------------------------------------

N = 2e2

P = 1e3

t0 = 0
t1 = 1

time_grid = seq( t0, t1, length.out = P )

Cov = exp_cov_function( time_grid, alpha = 0.3, beta = 0.4 )

Data_1 = generate_gauss_fdata( N, centerline = sin( 2 * pi * time_grid ), Cov = Cov )
Data_2 = generate_gauss_fdata( N, centerline = sin( 2 * pi * time_grid ), Cov = Cov )

mfD = mfData( time_grid, list( Data_1, Data_2 ) )

test_that( 'Kendall correlation with MEI ordering ',
           expect_silent( invisible( cor_spearman( mfD, ordering = 'MEI' ) ) ) )

test_that( 'Kendall correlation with MHI ordering',
           expect_silent( invisible( cor_spearman( mfD, ordering = 'MHI' ) ) ) )


# CASE STUDIES ------------------------------------------------------------

# ( With reference to the paper by Dalia Valencia, Rosa Lillo e Juan Romo)

N = 50
P = 50

t0 = 0
t1 = 1
time_grid = seq( t0, t1, length.out = P )

# Case 1
sigma_12 = 0.8

R = matrix( c( 1, sigma_12, sigma_12, 1 ), ncol = 2, nrow = 2 )

Z = matrix( rnorm( N * 2, 0, 1 ), ncol = 2, nrow = N ) %*% chol( R )

# X = t( ( t( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE )  ) + Z[ 1,] ) )^3
X = ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 1 ] )^3 +
    ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 1 ] )^2 +
    ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 1 ] ) * 3

Y = ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 2 ] )^2 +
    ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 2 ] ) * 7 / 8 +
    - 10

mfD = mfData( time_grid, list( X, Y ) )

cor_kendall( mfD, ordering = 'max' )
cor_kendall( mfD, ordering = 'area' )
cor_spearman( mfD, ordering = 'MEI' )
cor_spearman( mfD, ordering = 'MHI' )

# Case 2
sigma_12 = - 0.7

R = matrix( c( 1, sigma_12, sigma_12, 1 ), ncol = 2, nrow = 2 )

Z = matrix( rnorm( N * 2, 0, 1 ), ncol = 2, nrow = N ) %*% chol( R )

X = sin( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 1 ] )

Y = cos( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 2 ] )

mfD = mfData( time_grid, list( X, Y ) )

cor_kendall( mfD, ordering = 'max' )
cor_kendall( mfD, ordering = 'area' )
cor_spearman( mfD, ordering = 'MEI' )
cor_spearman( mfD, ordering = 'MHI' )


# Case 3
sigma_12 = 1

Z = rnorm( N, 0, 1)

X = ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z )^2

Y = ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z )^4

mfD = mfData( time_grid, list( X, Y ) )

test_that( ' Kendall correlation coeff with max ordering - Case 3',
           expect_equal( cor_kendall( mfD, ordering = 'max' ), 1 ) )
test_that( ' Kendall correlation coeff with area ordering - Case 3',
           expect_equal( cor_kendall( mfD, ordering = 'area' ), 1 ) )
test_that( ' Spearman correlation coeff with MEI ordering - Case 3',
           expect_equal( cor_spearman( mfD, ordering = 'MEI' ), 1 ) )
test_that( ' Spearman correlation coeff with MHI ordering - Case 3',
           expect_equal( cor_spearman( mfD, ordering = 'MHI' ), 1 ) )


# Case 4
sigma_12 = 1

Z = rnorm( N, 0, 1)

X = ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z )^2 +
  ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z ) * 7 +
  2

Y = ( ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z )^2 +
        ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z ) * 7 +
        2 )^3
mfD = mfData( time_grid, list( X, Y ) )

test_that( ' Kendall correlation coeff with max ordering - Case 4',
           expect_equal( cor_kendall( mfD, ordering = 'max' ), 1 ) )
test_that( ' Kendall correlation coeff with area ordering - Case 4',
           expect_equal( cor_kendall( mfD, ordering = 'area' ), 1 ) )
test_that( ' Spearman correlation coeff with MEI ordering - Case 4',
           expect_equal( cor_spearman( mfD, ordering = 'MEI' ), 1 ) )
test_that( ' Spearman correlation coeff with MHI ordering - Case 4',
           expect_equal( cor_spearman( mfD, ordering = 'MHI' ), 1 ) )

# Case 5
sigma_12 = 1

Z = rnorm( N, 0, 1)

X = ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z )^2 +
  ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z ) * 7 +
  2

Y = 1 - ( ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z )^2 +
        ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z ) * 7 +
        2 )^3
mfD = mfData( time_grid, list( X, Y ) )

test_that( ' Kendall correlation coeff with max ordering - Case 5',
           expect_equal( cor_kendall( mfD, ordering = 'max' ), -1 ) )
test_that( ' Kendall correlation coeff with area ordering - Case 5',
           expect_equal( cor_kendall( mfD, ordering = 'area' ), -1 ) )
test_that( ' Spearman correlation coeff with MEI ordering - Case 5',
           expect_equal( cor_spearman( mfD, ordering = 'MEI' ), -1 ) )
test_that( ' Spearman correlation coeff with MHI ordering - Case 5',
           expect_equal( cor_spearman( mfD, ordering = 'MHI' ), -1 ) )


# Case 6
sigma_12 = 0.6

R = matrix( c( 1, sigma_12, sigma_12, 1 ), ncol = 2, nrow = 2 )

Z = matrix( rnorm( N * 2, 0, 1 ), ncol = 2, nrow = N ) %*% chol( R )

X = exp( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 1 ] )

Y = ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 2 ] )^3 +
    ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 2 ] )^2 +
    ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 2 ] ) * 3

mfD = mfData( time_grid, list( X, Y ) )

cor_kendall( mfD, ordering = 'max' )
cor_kendall( mfD, ordering = 'area' )
cor_spearman( mfD, ordering = 'MEI' )
cor_spearman( mfD, ordering = 'MHI' )


# Case 7
sigma_12 = -0.8

R = matrix( c( 1, sigma_12, sigma_12, 1 ), ncol = 2, nrow = 2 )

Z = matrix( rnorm( N * 2, 0, 1 ), ncol = 2, nrow = N ) %*% chol( R )

X = exp( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 1 ] )^2

Y = cos( ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 2 ] ) )

mfD = mfData( time_grid, list( X, Y ) )

cor_kendall( mfD, ordering = 'max' )
cor_kendall( mfD, ordering = 'area' )
cor_spearman( mfD, ordering = 'MEI' )
cor_spearman( mfD, ordering = 'MHI' )


# Case 8
sigma_12 = 0.4

R = matrix( c( 1, sigma_12, sigma_12, 1 ), ncol = 2, nrow = 2 )

Z = matrix( rnorm( N * 2, 0, 1 ), ncol = 2, nrow = N ) %*% chol( R )

X = sin( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 1 ] )

Y = ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 2 ] )^2

mfD = mfData( time_grid, list( X, Y ) )

cor_kendall( mfD, ordering = 'max' )
cor_kendall( mfD, ordering = 'area' )
cor_spearman( mfD, ordering = 'MEI' )
cor_spearman( mfD, ordering = 'MHI' )

# Case 9
sigma_12 = 1

Z = rnorm( N, 0, 1 )

X = ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z )^2 +
    ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z ) * 9 +
    - 5

Y = cos( matrix( 3 * time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z )

mfD = mfData( time_grid, list( X, Y ) )

cor_kendall( mfD, ordering = 'max' )
cor_kendall( mfD, ordering = 'area' )
cor_spearman( mfD, ordering = 'MEI' )
cor_spearman( mfD, ordering = 'MHI' )

# Case 10
sigma_12 = 0.9

R = matrix( c( 1, sigma_12, sigma_12, 1 ), ncol = 2, nrow = 2 )

Z = matrix( rnorm( N * 2, 0, 1 ), ncol = 2, nrow = N ) %*% chol( R )

X = exp( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE )^2 + Z[ , 1 ] )

Y = ( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 2 ] )^2 +
      matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) * ( - 8 ) +
    ( matrix( 0, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 2 ] )

mfD = mfData( time_grid, list( X, Y ) )

cor_kendall( mfD, ordering = 'max' )
cor_kendall( mfD, ordering = 'area' )
cor_spearman( mfD, ordering = 'MEI' )
cor_spearman( mfD, ordering = 'MHI' )

# Case 11
sigma_12 = 0.

R = matrix( c( 1, sigma_12, sigma_12, 1 ), ncol = 2, nrow = 2 )

Z = matrix( rnorm( N * 2, 0, 1 ), ncol = 2, nrow = N ) %*% chol( R )

X = exp( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE )+ Z[ , 1 ] )

Y = sin( matrix( time_grid, nrow = N, ncol = P, byrow = TRUE ) + Z[ , 2 ] )

mfD = mfData( time_grid, list( X, Y ) )

cor_kendall( mfD, ordering = 'max' )
cor_kendall( mfD, ordering = 'area' )
cor_spearman( mfD, ordering = 'MEI' )
cor_spearman( mfD, ordering = 'MHI' )
