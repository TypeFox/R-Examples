library( "mvProbit" )
options( digits = 4 )

## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 10

# generate explanatory variables
xMat <- cbind( 
   const = rep( 1, nObs ),
   x1 = as.numeric( rnorm( nObs ) > 0 ),
   x2 = rnorm( nObs ),
   x3 = rnorm( nObs ) )

# model coefficients
beta <- cbind( c( -0.4,  0.8, -1.0, -1.4 ),
               c( -0.6, -0.5,  0.6, -1.2 ),
               c(  0.3,  0.9, -1.3,  0.8 ),
               c( -0.3, -1.0,  0.5, -0.8 ),
               c(  1.2, -0.6, -0.7,  1.1 ) )

# covariance matrix of error terms
sigma <- diag( 5 )
sigma[ upper.tri( sigma ) ] <- 
   c( 0.5, 0.4, -0.7, 0.8, -0.5, 0.8, 0, 0.6, -0.2, 0.6 )
sigma <- cov2cor( t( sigma ) %*% sigma )

# all parameters in a vector
allCoef <- c( c( beta ), sigma[ lower.tri( sigma ) ] )

# generate dependent variables
yMatLin <- xMat %*% beta 
yMat <- ( yMatLin + rmvnorm( nObs, sigma = sigma, pre0.9_9994 = TRUE ) ) > 0
colnames( yMat ) <- paste( "y", 1:5, sep = "" )
# (yMatLin > 0 )== yMat

# unconditional expectations of dependent variables
yExp <- mvProbitExp( ~ x1 + x2 + x3, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ) )
round( yExp, 3 )
yExpA <- mvProbitExp( ~ x1 + x2 + x3, coef = allCoef,
   data = as.data.frame( xMat ) )
all.equal( yExp, yExpA )
yExp2 <- pnorm( yMatLin )
all.equal( yExp, as.data.frame( yExp2 ) )


# conditional expectations of dependent variables
# (assuming that all other dependent variables are one)
yExpCond <- mvProbitExp( ~ x1 + x2 + x3, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = GenzBretz() )
round( yExpCond, 3 )
yExpCondA <- mvProbitExp( ~ x1 + x2 + x3, coef = allCoef,
   data = as.data.frame( xMat ), cond = TRUE, algorithm = GenzBretz() )
all.equal( yExpCond, yExpCondA )
yExpCond2 <- matrix( NA, nrow = nObs, ncol = ncol( yMat ) )
for( i in 1:nObs ) {
   for( k in 1:ncol( yMat ) ) {
      set.seed( 123 )
      numerator <- pmvnorm( upper = yMatLin[ i, ], sigma = sigma )
      set.seed( 123 )
      denominator <- pmvnorm( upper = yMatLin[ i, -k ], sigma = sigma[ -k, -k ] )
      yExpCond2[ i, k ] <- numerator / denominator
   }
}
all.equal( yExpCond, as.data.frame( yExpCond2 ) )


# conditional expectations of dependent variables
# (assuming that all other dependent variables are as observed)
yExpCondObs <- mvProbitExp( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ), 
   cond = TRUE, algorithm = GenzBretz() )
round( yExpCondObs, 3 )
yExpCondObsA <- mvProbitExp( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3, 
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) ), 
   cond = TRUE, algorithm = GenzBretz() )
all.equal( yExpCond, yExpCondA )
yExpCondObs2 <- matrix( NA, nrow = nObs, ncol = ncol( yMat ) )
for( i in 1:nObs ){
   for( k in 1:ncol( yMat ) ) {
      ySign <- 2 * yMat[ i, ] - 1
      ySign[ k ] <- 1
      yLinTmp <- yMatLin[ i, ] * ySign
      sigmaTmp <- diag( ySign ) %*% sigma %*% diag( ySign )
      set.seed( 123 )
      numerator <- pmvnorm( upper = yLinTmp, sigma = sigmaTmp )
      set.seed( 123 )
      denominator <- pmvnorm( upper = yLinTmp[ -k ], sigma = sigmaTmp[ -k, -k ] )
      yExpCondObs2[ i, k ] <- numerator / denominator
   }
}
colnames( yExpCondObs2 ) <- names( yExpCondObs )
all.equal( yExpCondObs, as.data.frame( yExpCondObs2 ) )


# calculating log likelihood value(s)
logLikVal <- mvProbitLogLik( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   algorithm = GenzBretz() )
round( logLikVal, 3 )
logLikValA <- mvProbitLogLik( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3, 
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) ),
   algorithm = GenzBretz() )
all.equal( logLikVal, logLikValA )

# calculating log likelihood value(s) with one-sided gradients
logLikValGrad1 <- mvProbitLogLik( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3,
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   oneSidedGrad = TRUE, algorithm = GenzBretz() )
round( c( logLikValGrad1 ), 3 )
round( attr( logLikValGrad1, "gradient" ), 3 )
logLikValGrad1A <- mvProbitLogLik( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3,
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) ),
   oneSidedGrad = TRUE, algorithm = GenzBretz() )
all.equal( logLikValGrad1, logLikValGrad1A )

# calculating log likelihood value(s) with two-sided gradients
logLikValGrad <- mvProbitLogLik( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3,
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   returnGrad = TRUE, algorithm = GenzBretz() )
round( c( logLikValGrad ), 3 )
round( attr( logLikValGrad, "gradient" ), 3 )
# now manually
llTmp <- function( coef ) {
   betaTmp <- coef[ 1:20 ]
   sigmaTmp <- diag( 5 )
   sigmaTmp[ lower.tri( sigmaTmp ) ] <- coef[ -(1:20) ]
   sigmaTmp[ upper.tri( sigmaTmp ) ] <- t( sigmaTmp )[ upper.tri( sigmaTmp ) ]
   result <- mvProbitLogLik( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3,
      coef = betaTmp, sigma = sigmaTmp, 
      data = as.data.frame( cbind( xMat, yMat ) ), algorithm = GenzBretz() )
   return( result )
}
logLikValGrad2 <- numericGradient( llTmp, allCoef )
round( logLikValGrad2, 3 )
# attr( logLikValGrad1, "gradient" ) / logLikValGrad2 - 1
# attr( logLikValGrad1, "gradient" ) - logLikValGrad2
all.equal( attr( logLikValGrad1, "gradient" ), logLikValGrad2,
  tol = 1e-5, check.attributes = FALSE )
# attr( logLikValGrad, "gradient" ) / logLikValGrad2 - 1
# attr( logLikValGrad, "gradient" ) - logLikValGrad2
all.equal( attr( logLikValGrad, "gradient" ), logLikValGrad2,
  check.attributes = FALSE )


# calculating marginal effects, unconditional
margEffUnc <- mvProbitMargEff( ~ x1 + x2 + x3, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), vcov = diag( 30 ),
   returnJacobian = TRUE )
round( margEffUnc, 3 )
round( attr( margEffUnc, "vcov" )[ 1:3, , ], 2 )
round( drop( attr( margEffUnc, "vcov" )[ nObs, , ] ), 2 )
round( attr( margEffUnc, "jacobian" )[ 1:3, , ], 2 )
round( drop( attr( margEffUnc, "jacobian" )[ nObs, , ] ), 2 )
print( summary( margEffUnc ), digits = c( 3, 3, 2, 2, 2 ) )
margEffUncA <- mvProbitMargEff( ~ x1 + x2 + x3, coef = allCoef,
   data = as.data.frame( xMat ), vcov = diag( 30 ), returnJacobian = TRUE )
all.equal( margEffUnc, margEffUncA )

# calculating marginal effects, conditional
# (assuming that all other dependent variables are one)
margEffCond <- mvProbitMargEff( ~ x1 + x2 + x3, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = GenzBretz() )
round( margEffCond, 3 )
margEffCondA <- mvProbitMargEff( ~ x1 + x2 + x3, coef = allCoef,
   data = as.data.frame( xMat ), cond = TRUE,
   algorithm = GenzBretz() )
all.equal( margEffCond, margEffCondA )
# now with variance covariance matrix
margEffCondV <- mvProbitMargEff( ~ x1 + x2 + x3, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat )[ c(1,5,10), ], cond = TRUE, 
   vcov = diag( 30 ), returnJacobian = TRUE, algorithm = GenzBretz() )
round( attr( margEffCondV, "vcov" ), 2 )
round( drop( attr( margEffCondV, "vcov" )[ 1, , ] ), 2 )
round( attr( margEffCondV, "jacobian" ), 2 )
round( drop( attr( margEffCondV, "jacobian" )[ 1, , ] ), 2 )
print( summary( margEffCondV ), digits = c( 3, 3, 2, 2, 2 ) )

# calculating marginal effects, conditional
# (assuming that all other dependent variables are as observed)
margEffCondObs <- mvProbitMargEff( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3,
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ), 
   cond = TRUE, algorithm = GenzBretz() )
round( margEffCondObs, 3 )
margEffCondObsA <- mvProbitMargEff( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3,
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) ), cond = TRUE,
   algorithm = GenzBretz() )
all.equal( margEffCondObs, margEffCondObsA )
# now with variance covariance matrix
margEffCondObsV <- mvProbitMargEff( cbind( y1, y2, y3, y4, y5 ) ~ x1 + x2 + x3, 
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) )[ c(1,5,10), ], 
   cond = TRUE, vcov = diag( 30 ), returnJacobian = TRUE, algorithm = GenzBretz() )
round( attr( margEffCondObsV, "vcov" ), 2 )
round( drop( attr( margEffCondObsV, "vcov" )[ 1, , ] ), 2 )
round( attr( margEffCondObsV, "jacobian" ), 2 )
round( drop( attr( margEffCondObsV, "jacobian" )[ 1, , ] ), 2 )
print( summary( margEffCondObsV ), digits = c( 3, 3, 2, 2, 2 ) )

