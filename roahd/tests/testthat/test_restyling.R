


# RESTYLING UNIVARIATE DEPTHS ---------------------------------------------

###### BD

N = 1e2

P = 1e2


grid = seq( 0, 1, length.out = 1e2 )

centerline = sin( 2 * pi * grid )
Cov = exp_cov_function( grid, alpha = 0.2, beta = 0.3 )

Data = generate_gauss_fdata( N,
                             centerline,
                             Cov )
fD = fData( grid, Data )

test_that( 'Restyling test on BD ',
           expect_identical( BD( fD ),
                             BD( Data ) ) )

###### MBD

N = 1e2

P = 1e2


grid = seq( 0, 1, length.out = 1e2 )

centerline = sin( 2 * pi * grid )
Cov = exp_cov_function( grid, alpha = 0.2, beta = 0.3 )

Data = generate_gauss_fdata( N,
                             centerline,
                             Cov )
fD = fData( grid, Data )

test_that( 'Restyling test on MBD - 1 ',
           expect_identical( MBD( fD, manage_ties = TRUE ),
                             MBD( Data, manage_ties = TRUE ) ) )

test_that( 'Restyling test on MBD - 2 ',
           expect_identical( MBD( fD, manage_ties = FALSE ),
                             MBD( Data, manage_ties = FALSE ) ) )

###### RELATIVE BD and RELATIVE MBD

N = 1e2

P = 1e2

grid = seq( 0, 1, length.out = 1e2 )

centerline_1 = sin( 2 * pi * grid )
centerline_2 = sin( 2 * pi * grid ) + 0.5

Cov = exp_cov_function( grid, alpha = 0.2, beta = 0.3 )

Data_reference = generate_gauss_fdata( N,
                             centerline_1,
                             Cov )
Data_target = generate_gauss_fdata( N,
                                    centerline_2,
                                    Cov )
fD_target = fData( grid, Data_target )
fD_reference = fData( grid, Data_reference )

test_that( 'Restyling test on BD ',
           expect_identical( BD_relative( fD_target, fD_reference ),
                             BD_relative( Data_target, Data_reference ) ) )
test_that( 'Restyling test on MBD ',
           expect_identical( MBD_relative( fD_target, fD_reference ),
                             MBD_relative( Data_target, Data_reference ) ) )


# RESTYLING INDEXES -------------------------------------------------------

#### EI and MEI
N = 1e2

P = 1e2

grid = seq( 0, 1, length.out = 1e2 )

centerline = sin( 2 * pi * grid )

Cov = exp_cov_function( grid, alpha = 0.2, beta = 0.3 )

Data = generate_gauss_fdata( N,
                             centerline,
                             Cov )
fD = fData( grid, Data )

test_that( 'Restyling test on EI',
           expect_identical( EI( fD ),
                             EI( Data ) ) )

test_that( 'Restyling test on MEI',
           expect_identical( MEI( fD ),
                             MEI( Data ) ) )

#### HI and MHI
test_that( 'Restyling test on EI',
           expect_identical( HI( fD ),
                             HI( Data ) ) )

test_that( 'Restyling test on MHI',
           expect_identical( MHI( fD ),
                             MHI( Data ) ) )

#### HRD and MHRD
test_that( 'Restyling test on EI',
           expect_identical( HRD( fD ),
                             HRD( Data ) ) )

test_that( 'Restyling test on MEI',
           expect_identical( MHRD( fD ),
                             MHRD( Data ) ) )
