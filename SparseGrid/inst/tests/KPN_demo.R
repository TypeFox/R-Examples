# Copyright (C) 2011 Jelmer Ypma. All Rights Reserved.
# This code is published under the GPL.
#
# File:   KPN_demo.R
# Author: Jelmer Ypma
# Date:   21 April 2012
#
# Approximate second moments for multivariate normal distribution, by integrating
#    \int_{-\infty}^{\infty} (Cx)(Cx)' f(x) dx.
# Where f(x) is the multivariate normal distribution, with mu = 0, and Sigma = I.
# g(x) = (Cx)(Cx)', is the function that want to integrate over.
#
# If x is normal with variance-covariance matrix I, then Cx is normal with 
# variance-covariance matrix CC', which means we can easily calculate the true
# value in this example, since we're integrating over the second moments, which
# should give the variance-covariance matrix by definition.

library('SparseGrid')

# accuracy
k <- 2

# Cholesky decomposition of variance-covariance matrix  
mat_C  <- rbind( c( 1, 0, 0 ), c( .5, 1.5, 0 ), c( .3, .7, 2 ) )

# Determine dimension of C, and calculate Sigma = CC'
dimension <- nrow( mat_C )
varcov.true <- mat_C %*% t( mat_C )      # variance-covariance matrix

# define integration grid
tmp.grid <- createSparseGrid('KPN', dimension=dimension, k=k)

# integrand: some function that evaluates g(x): (R times D)->(R times 1)
func <- function( x, mat_C ) {
    y <- mat_C %*% x        # create correlated variables
    return( y %*% t(y) )    # outer product (y is column vector)
}

# approximate integral
varcov.approx <- matrix( apply( tmp.grid$nodes, 1, func, mat_C=mat_C ) %*% tmp.grid$weights, nrow=dimension, ncol=dimension )

varcov.true
varcov.approx
varcov.true - varcov.approx
