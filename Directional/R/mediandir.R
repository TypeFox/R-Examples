################################
#### Median direction
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
#### References: Fisher, N. I. (1985). Spherical medians.
#### Journal of the Royal Statistical Society. Series B, 47(2): 342-348.
#### Fisher, N. I., Lewis, T., & Embleton, B. J. (1987).
#### Statistical analysis of spherical data. Cambridge university press.
################################

mediandir = function(x) {
  ## x is the directional data

  x = as.matrix(x)
  x = x / sqrt( rowSums(x^2) )
  n = nrow(x)  ;  p = ncol(x)
  funa = function(pa) {
    pa = pa / sqrt( sum(pa^2) )
    f = mean( acos( x %*% t( t(pa) ) ) )
    f
  }

  bar = optim( colMeans(x), funa, control = list(maxit = 10000) )
  bar = optim( bar$par, funa, control = list(maxit = 10000) )
  bar = optim( bar$par, funa, control = list(maxit = 10000) )
  bar = optim( bar$par, funa, control = list(maxit = 10000) )
  med = bar$par
  med / sqrt( sum(med^2) )

}
