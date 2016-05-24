
get.x11 <- function( cv, x, i, j ) x[ cv$row != i, cv$col != j ]
get.x12 <- function( cv, x, i, j ) x[ cv$row != i, cv$col == j ]
get.x21 <- function( cv, x, i, j ) x[ cv$row == i, cv$col != j ]
get.x22 <- function( cv, x, i, j ) x[ cv$row == i, cv$col == j ]

check <- function( cv, x )
{
    expected <- matrix( 0, cv$krow * cv$kcol, cv$maxrank + 1 )
    for( i in 1:cv$krow )
    {
        for( j in 1:cv$kcol )
        {
            x11 <- get.x11( cv, x, i, j );
            x11.svd <- svd( x11, nu=cv$maxrank, nv=cv$maxrank )
            u1 <- x11.svd$u
            v1 <- x11.svd$v
            d  <- x11.svd$d
            
            x12 <- get.x12( cv, x, i, j )
            x21 <- get.x21( cv, x, i, j )
            x22 <- get.x22( cv, x, i, j )
            
            row <- i + ( j-1 )*cv$krow
            expected[ row, 1 ] <- mean( x22^2 )
            for( k in seq_len( cv$maxrank ) )
            {
                x22.hat <- (     ( x21 %*% v1[,1:k,drop=FALSE] )
                             %*% diag( 1.0 / d[1:k], k, k )
                             %*% ( t( u1[,1:k,drop=FALSE] ) %*% x12 ) )
                expected[ row, k + 1 ] <- mean( ( x22 - x22.hat )^2 )                
            }
        }
    }
    
    actual <- cv$msep
    
    if( !all( abs( expected - actual )/expected < 1e-8 ) )
    {
        cat( "Expected:\n")
        print( expected )
        
        cat( "Actual:\n")
        print( actual )
        TRUE
    }
    else
        FALSE
}

rmatrix <- function( size )
{
    M <- floor( runif( 1, 4, max( 5, size ) ) )
    N <- floor( runif( 1, 4, max( 5, size ) ) )
    x <- matrix( rnorm( M*N ), M, N )
    x
}

rfold <- function( n )
{
    floor( runif( 1, 2, max( 3, n ) ) )
}

rmaxrank <- function( n, p, krow, kcol )
{
    maxrank <- min( n*( 1 - 1/krow ) , p*( 1 - 1/kcol ) )
    floor( runif( 1, 0, maxrank + 1 ) )
}

rtest <- function( size )
{
    x       <- rmatrix( size )
    krow    <- rfold( nrow( x ) )
    kcol    <- rfold( ncol( x ) )
    maxrank <- rmaxrank(  nrow( x ), ncol( x ), krow, kcol )

    list( x=x, krow=krow, kcol=kcol, maxrank=maxrank )
}

sizes <- function( ntests, each=2 )
{
    if( ntests > 0 )
        rep(4 + sqrt( seq( 0, ntests/each, length=ceiling( ntests/each ) ) )
                    , each=each )[ 1:ntests ]
    else
        c()
}

require( bcv )
ntests <- 500
s <- sizes( ntests )
set.seed( 0 )

nsuc <- 0
for (size in s)
{
    cat( '.' )
    test <- rtest( size )
    
    cv <- suppressWarnings( 
              cv.svd.gabriel( test$x, test$krow, test$kcol, test$maxrank ) )
    if( !check( cv, test$x ) )
        nsuc <- nsuc + 1
}

if( nsuc == ntests ) {
    cat("Passed", nsuc, "tests!\n")
} else {
    cat("Not all tests passed.\n")
}

