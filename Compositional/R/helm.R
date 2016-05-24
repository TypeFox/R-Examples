################################
#### The Helmert sub-matrix
#### Tsagris Michail 5/2011  
#### References: John Aitchison (2003) 
#### The Statistical Analysis of Compositional Data p. 99 Blackburn Press 
#### Lancaster H. O. (1965). The Helmert matrices. 
#### The American Mathematical Monthly 72(1): 4-12.
################################

helm <- function(n) {
  mat <- matrix(0, n - 1, n) 
  i <- 2:n
  r <- 1 / sqrt( i * (i - 1) )
  for ( j in 1:(n - 1 ) ) { 
    mat[j, 1: c(j + 1) ] <- c( rep(r[j], j), -( j * r[j] ) )
  }
  mat
}