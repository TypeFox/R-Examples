## -----------------------------------------------------------------------------
## Fonction runifSphere
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

runifSphere = function(dimension,N,radius) {
    tmp = matrix(rnorm(dimension*N), dimension, N)
    tmp <- t(tmp)*c(1/sqrt(rep(1,dim(tmp)[1])%*%tmp^2))*runif(N)^(1/dimension)*radius
    colnames(tmp) <- rep(c('x', 'y'), length.out = dimension)
    return(tmp)
}