################################
#### Median direction
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
#### References: Fisher, N. I. (1985). Spherical medians.
#### Journal of the Royal Statistical Society. Series B, 47(2): 342-348.
#### Fisher, N. I., Lewis, T., & Embleton, B. J. (1987).
#### Statistical analysis of spherical data. Cambridge university press.
#### Cabrera, J., & Watson, G. S. (1990). On a spherical median related distribution.
#### Communications in Statistics-Theory and Methods, 19(6): 1973-1986.
################################

mediandir_2 = function(x) {
  ## x is the directional data

  x = as.matrix(x)
  x = x / sqrt( rowSums(x^2) )
  n = nrow(x)  ;  p = ncol(x)
  pa = colMeans(x)
  u <- matrix(nrow = 1000, ncol = p)
  u[1, ] <- pa / sqrt( sum(pa^2) )
  ww <- as.vector( sqrt( 1 - ( x %*% u[1, ] )^2 ) )
  u[2, ] <- colSums( x / ww )
  u[2, ] <- u[2, ] / sqrt( sum(u[2, ]^2) )
  i <- 2

  while ( sum( abs(u[i, ] - u[i - 1, ]) ) > 1e-10 ) {
    i <- i +1
    ww <- as.vector( sqrt( 1 - ( x %*% u[i - 1, ] )^2 ) )
    u[i, ] <- colSums(x / ww )
    u[i, ] <- u[i, ] / sqrt( sum(u[i, ]^2) )
  }

  u[i, ]
}
