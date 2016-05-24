GaussianMatrix <- function(N, M) {
 # Part of R1Magic by mehmet.suzen@physics.org
 return( matrix( rnorm(N * M ),  M,N ) )
}

