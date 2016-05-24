##
##  m e s h g r i d . R Test suite
##


meshgrid <- pracma::meshgrid

identical(meshgrid(1:3, 10:14)$X,
          matrix(rep(c(1:3), each = 5), nr = 5, nc = 3))
identical(meshgrid(1:3, 10:14)$Y, 
          matrix(rep(10:14, times = 3), nr = 5, nc = 3))
