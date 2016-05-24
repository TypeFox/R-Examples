

# TESTING HI AND MHI ------------------------------------------------------

N = 20

time_grid = seq( 0, 1, length.out = N * 1e2 )

Data = matrix( 0, nrow = N, ncol = length( time_grid ) )

for( iObs in 1 : N )
{
  Data[ iObs, ] =  as.numeric( time_grid >= ( iObs - 1 ) / N & time_grid < iObs / N )
}

Data[ N, length( time_grid ) ] = 1

test_that( "Correctness of HI",
           expect_equal( HI( Data ), rep( ( N - 1 ) / N, N ) ) )

test_that( "Correctness of MHI",
           expect_equal( MHI( Data ), rep( ( N + ( N - 1 )^2 ) / N^2 , N ) ) )
