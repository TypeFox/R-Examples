orthonormal.polynomials <- function( recurrences, p.0 )
{
###
### This function returns a list with n+1 elements
### containing the order orthonormal polynomials for orders k=0,1,...,n
###
### Parameter
### recurrences = a data frame containing the parameters c, d, e and f
### p.0 = a polynomial object for the order 0 orthonormal polynomial
###
    require( polynom )
    np1 <- nrow( recurrences )
    n <- np1 - 1
    c <- recurrences$c
    d <- recurrences$d
    e <- recurrences$e
    f <- recurrences$f
    polynomials <- as.list( rep( NULL, np1 ) )
    polynomials[[1]] <- p.0
    j <- 0
    while ( j < n ) {
        cj <- c[j+1]
        dj <- d[j+1]
        ej <- e[j+1]
        fj <- f[j+1]
        monomial <- polynomial( c( dj, ej ) )
        if ( j == 0 ) {
            p.jp1 <- ( monomial * p.0 ) / cj
        }
        else {
            p.jm1 <- polynomials[[j]]
            p.j   <- polynomials[[j+1]]
            p.jp1 <- ( monomial * p.j - fj * p.jm1 ) / cj
        }
        polynomials[[j+2]] <- p.jp1
        j <- j + 1
    }
    return( polynomials )
}
