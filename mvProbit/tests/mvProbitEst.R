library( "mvProbit" )
options( digits = 4 )

## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 50

# generate explanatory variables
xMat <- cbind( 
   const = rep( 1, nObs ),
   x1 = as.numeric( rnorm( nObs ) > 0 ),
   x2 = rnorm( nObs ) )

# model coefficients
beta <- cbind( c(  0.8,  1.2, -0.8 ),
               c( -0.6,  1.0, -1.6 ),
               c(  0.5, -0.6,  1.2 ) )

# covariance matrix of error terms
sigma <- miscTools::symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )

# generate dependent variables
yMatLin <- xMat %*% beta 
yMat <- ( yMatLin + rmvnorm( nObs, sigma = sigma, pre0.9_9994 = TRUE ) ) > 0
colnames( yMat ) <- paste( "y", 1:3, sep = "" )

# create data frame
dat <- as.data.frame( cbind( xMat, yMat ) )

# estimation with the BHHH algorithm, two-sided gradients
estResultBHHH <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma,
   data = dat, tol = 0.5,
   algorithm = GenzBretz() )
print( estResultBHHH )
summary( estResultBHHH )
logLik( estResultBHHH )
estResultBHHHA <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta, sigma[ lower.tri( sigma ) ] ),
   data = dat, tol = 0.5,
   algorithm = GenzBretz() )
all.equal( estResultBHHH, estResultBHHHA )

# estimation with the BHHH algorithm, one-sided gradients
estResultBHHH1 <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma,
   data = dat, tol = 0.5,
   algorithm = GenzBretz(), oneSidedGrad = TRUE )
print( estResultBHHH1 )
summary( estResultBHHH1 )
logLik( estResultBHHH1 )
all.equal( estResultBHHH, estResultBHHH1, tol = 1e-5 )

# estimation with the BFGS algorithm, two-sided gradients
estResultBFGS <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, method = "BFGS",
   data = dat, 
   reltol = 0.5, algorithm = GenzBretz() )
print( estResultBFGS )
summary( estResultBFGS )
logLik( estResultBFGS )

# estimation with the BFGS algorithm, one-sided gradients
estResultBFGS1 <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, method = "BFGS", 
   data = dat, 
   reltol = 0.5, algorithm = GenzBretz(), oneSidedGrad = TRUE )
print( estResultBFGS1 )
summary( estResultBFGS1 )
logLik( estResultBFGS1 )
all.equal( estResultBFGS, estResultBFGS1, tol = 1e-5 )

# estimation with the BFGS algorithm, one-sided gradients, no starting values
estResultBFGS1a <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   data = dat, method = "BFGS",
   reltol = 0.5, algorithm = GenzBretz(), oneSidedGrad = TRUE )
print( estResultBFGS1a )
summary( estResultBFGS1a )
logLik( estResultBFGS1a )

# estimation with the BFGS algorithm, one-sided gradients, no starting values for beta
estResultBFGS1b <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   startSigma = sigma, data = dat, 
   method = "BFGS", reltol = 0.5, algorithm = GenzBretz(), oneSidedGrad = TRUE )
print( estResultBFGS1b )
summary( estResultBFGS1b )
logLik( estResultBFGS1b )

# estimation with the BFGS algorithm, one-sided gradients, no starting values for sigma
estResultBFGS1s <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), data = dat, 
   method = "BFGS", reltol = 0.5, algorithm = GenzBretz(), oneSidedGrad = TRUE )
print( estResultBFGS1s )
summary( estResultBFGS1s )
logLik( estResultBFGS1s )

# estimation with the BFGS algorithm, Miwa algorithm for obtaining integrals
estResultBFGSm <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, 
   data = dat, method = "BFGS",
   reltol = 0.5, algorithm = Miwa( steps = 64 ) )
print( estResultBFGSm )
summary( estResultBFGSm )
logLik( estResultBFGSm )
all.equal( estResultBFGS, estResultBFGSm, tol = 1e-3 )

# estimation with the BFGS algorithm, GHK algorithm for obtaining integrals
estResultBFGSg <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, 
   data = dat, method = "BFGS",
   reltol = 0.5 )
print( estResultBFGSg )
summary( estResultBFGSg )
logLik( estResultBFGSg )
all.equal( estResultBFGS, estResultBFGSg, tol = 1e-2 )
all.equal( estResultBFGSm, estResultBFGSg, tol = 1e-2 )

# estimation with the Nelder-Mead algorithm
estResultNM <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, 
   data = dat, 
   method = "NM", reltol = 0.05, algorithm = GenzBretz() )
print( estResultNM )
summary( estResultNM )
logLik( estResultNM )

# estimation with NA in an explanatory variable
dat$x1Na <- dat$x1
dat$x1Na[7] <- NA
estResultNax <- try( mvProbit( cbind( y1, y2, y3 ) ~ x1Na + x2,
   data = dat, method = "BFGS", reltol = 0.5, algorithm = GenzBretz() ) )
estResultNaxStart <- try( mvProbit( cbind( y1, y2, y3 ) ~ x1Na + x2,
   start = c( beta ), startSigma = sigma,
   data = dat, method = "BFGS", reltol = 0.5, algorithm = GenzBretz() ) )

# estimation with NA in a dependant variable
dat$y2Na <- dat$y2
dat$y2Na[9] <- NA
estResultNay <- try( mvProbit( cbind( y1, y2Na, y3 ) ~ x1 + x2,
   data = dat, method = "BFGS", reltol = 0.5, algorithm = GenzBretz() ) )

# estimation with NA both in a dependant variable and in an explanatory variable
estResultNaxy <- try( mvProbit( cbind( y1, y2Na, y3 ) ~ x1Na + x2,
   start = c( beta ), startSigma = sigma,
   data = dat, method = "BFGS", reltol = 0.5, algorithm = GenzBretz() ) )

# estimation with infinity in an explanatory variable
dat$x2Inf <- dat$x2
dat$x2Inf[15] <- Inf
estResultInf <- try( mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2Inf,
   data = dat, method = "BFGS", reltol = 0.5, algorithm = GenzBretz() ) )
estResultInfStart <- try( mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2Inf,
   start = c( beta ), startSigma = sigma,
   data = dat, method = "BFGS", reltol = 0.5, algorithm = GenzBretz() ) )

# estimation with a factor as explanatory variable (x1 with 2 levels)
dat$x1Fac <- as.factor( ifelse( dat$x1 == 0, "green", "red" ) )
estResultFac <- mvProbit( cbind( y1, y2, y3 ) ~ x1Fac + x2,
   start = c( beta ), startSigma = sigma,
   data = dat, method = "BFGS", reltol = 0.5, algorithm = GenzBretz() )
print( estResultFac )
summary( estResultFac )
logLik( estResultFac )
all.equal( estResultBFGS, estResultFac )

# estimation with a factor as explanatory variable (x1 with 3 levels)
dat$x1Fac3 <- as.factor(
   ifelse( rnorm( nObs ) <= -0.5, "brown", as.character( dat$x1Fac ) ) )
estResultFac13 <- mvProbit( cbind( y1, y2, y3 ) ~ x1Fac3 + x2,
   data = dat, method = "BFGS", reltol = 0.5, algorithm = GenzBretz() )
print( estResultFac13 )
summary( estResultFac13 )
logLik( estResultFac13 )

# estimation with a factor as explanatory variable (x2 with 3 levels)
dat$x2Fac3 <- as.factor( ifelse( dat$x2 < -0.5, "low", 
   ifelse( dat$x2 > 0.5, "high", "medium" ) ) )
estResultFac23 <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2Fac3,
   data = dat, method = "BFGS", reltol = 0.5, algorithm = GenzBretz() )
print( estResultFac23 )
summary( estResultFac23 )
logLik( estResultFac23 )

# estimation with a factor as explanatory variable (x3 with 3 levels)
dat$x3Fac3 <- as.factor( ifelse( rnorm( nObs ) < -0.5, "low", 
   ifelse( rnorm( nObs ) < 0, "medium", "high" ) ) )
estResultFac33 <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2 + x3Fac3,
   data = dat, method = "BFGS", reltol = 0.5, algorithm = GenzBretz() )
print( estResultFac33 )
summary( estResultFac33 )
logLik( estResultFac33 )


## testing the logLik method
# argument 'coef'
all.equal( logLik( estResultBHHH ), 
   logLik( estResultBHHH, coef = coef( estResultBHHH ) ) )
logLik( estResultBHHH ) -
   logLik( estResultBHHH, coef = coef( estResultBHHH ) * 0.999 )

# argument 'data'
all.equal( logLik( estResultBHHH ), 
   logLik( estResultBHHH, data = dat ) )
logLik( estResultBHHH ) -
   logLik( estResultBHHH, data = as.data.frame( cbind( xMat * 0.999, yMat ) ) )

# argument 'algorithm'
all.equal( logLik( estResultBHHH ), 
   logLik( estResultBHHH, algorithm = GenzBretz() ) )
logLik( estResultBHHH ) -
   logLik( estResultBHHH, algorithm = Miwa() )

# argument 'nGHK'
all.equal( logLik( estResultBHHH ), 
   logLik( estResultBHHH, nGHK = 5555 ) )
logLik( estResultBHHH, algorithm = "GHK" ) -
   logLik( estResultBHHH, algorithm = "GHK", nGHK = 2000 )

# argument 'random.seed'
all.equal( logLik( estResultBHHH ), 
   logLik( estResultBHHH, random.seed = 123 ) )
logLik( estResultBHHH ) -
   logLik( estResultBHHH, random.seed = 1234 )

# restore original data frame without factors, 
# because otherwise some of the following commands fail, 
# as the mean of factor variables cannot be calculated
datFac <- dat
dat <- as.data.frame( cbind( xMat, yMat ) )

# marginal effects based on estimated coefficients with covariance matrix
# unconditional marginal effects (with Jacobian)
margEffUnc <- margEff( estResultBFGS, calcVCov = TRUE, returnJacobian = TRUE )
round( margEffUnc, 3 )
round( attr( margEffUnc, "vcov" )[ 1:5, , ], 2 )
round( drop( attr( margEffUnc, "vcov" )[ nObs, , ] ), 2 )
round( attr( margEffUnc, "jacobian" )[ 1:5, , ], 2 )
round( drop( attr( margEffUnc, "jacobian" )[ nObs, , ] ), 2 )
print( summary( margEffUnc ), digits = c( 3, 3, 2, 2, 2 ) )
# now with explicitly specifying dummy variables
margEffUncD <- margEff( estResultBFGS, calcVCov = TRUE, returnJacobian = TRUE,
   dummyVars = c( "x1" ) )
all.equal( margEffUncD, margEffUnc )
# now with seemingly no dummy variables
margEffUncD0 <- margEff( estResultBFGS, calcVCov = TRUE, returnJacobian = TRUE,
   dummyVars = NULL )
print( summary( margEffUncD0 ), digits = c( 3, 3, 2, 2, 2 ) )
# now with seemingly only dummy variables
margEffUncDA <- margEff( estResultBFGS, calcVCov = TRUE, returnJacobian = TRUE,
   dummyVars = c( "x1", "x2" ) )
print( summary( margEffUncDA ), digits = 3 )
# now with returned Jacobian but without variance covariance matrix
margEffUncJac <- margEff( estResultBFGS, returnJacobian = TRUE )
all.equal( attr( margEffUncJac, "jacobian" ), attr( margEffUnc, "jacobian" ) )
# now including mean values of the marginal effects
margEffUncM <- margEff( estResultBFGS, calcVCov = TRUE, returnJacobian = TRUE,
   addMean = TRUE )
all.equal( margEffUnc, margEffUncM[ -(nObs+1), ], check.attributes = FALSE )
round( margEffUncM[ nObs:(nObs+1), ], 3 )
all.equal( attr( margEffUnc, "vcov" ), attr( margEffUncM, "vcov" )[ 1:nObs, , ] )
round( attr( margEffUncM, "vcov" )[ nObs:(nObs+1), , ], 2 )
round( drop( attr( margEffUncM, "vcov" )[ nObs+1, , ] ), 2 )
all.equal( attr( margEffUnc, "jacobian" ), attr( margEffUncM, "jacobian" )[ 1:nObs, , ] )
round( attr( margEffUncM, "jacobian" )[ nObs:(nObs+1), , ], 2 )
round( drop( attr( margEffUncM, "jacobian" )[ nObs+1, , ] ), 2 )
all.equal( summary( margEffUnc )[ , ], 
   summary( margEffUncM )[ 1:( 6 * nObs ), ], check.attributes = FALSE )
printCoefmat( round(
  summary( margEffUncM )[ -( 1:( 6 * ( nObs - 1 ) ) ), ], digits = 3 ),
  digits = 3 )
# now at mean values of explanatory variables
margEffUncMean <- margEff( estResultBFGS, calcVCov = TRUE, 
   data = as.data.frame( t( colMeans( xMat ) ) ) )
print( summary( margEffUncMean ), digits = c( 3, 3, 2, 2, 2 ) )
# now with argument 'atMean'
margEffUncMeanA <- margEff( estResultBFGS, calcVCov = TRUE, atMean = TRUE )
all.equal( margEffUncMeanA, margEffUncMean )

# conditional marginal effects
# (assuming that all other dependent variables are as observed)
margEffCondObs <- margEff( estResultBFGS, cond = TRUE, algorithm = GenzBretz() )
round( margEffCondObs, 3 )

# conditional marginal effects with covariance matrix at sample mean
# (assuming that all other dependent variables are at there modal values)
# (with Jacobian)
margEffCondObsCov <- margEff( estResultBFGS, cond = TRUE,
   atMean = TRUE, othDepVar = c( colMedians( yMat * 1 ) ), 
   algorithm = GenzBretz(), calcVCov = TRUE, returnJacobian = TRUE )
round( margEffCondObsCov, 3 )
round( attr( margEffCondObsCov, "vcov" ), 2 )
round( drop( attr( margEffCondObsCov, "vcov" ) ), 2 )
round( attr( margEffCondObsCov, "jacobian" ), 2 )
round( drop( attr( margEffCondObsCov, "jacobian" ) ), 2 )
print( summary( margEffCondObsCov ), digits = c( 3, 3, 2, 2, 2 ) )
# now with explicitly specifying dummy variables
margEffCondObsCovD <- margEff( estResultBFGS, cond = TRUE,
   atMean = TRUE, othDepVar = c( colMedians( yMat * 1 ) ),
   algorithm = GenzBretz(),
   calcVCov = TRUE, returnJacobian = TRUE, dummyVars = c( "x1" ) )
all.equal( margEffCondObsCovD, margEffCondObsCov )
# now with seemingly no dummy variables
margEffCondObsCovD0 <- margEff( estResultBFGS, cond = TRUE,
   atMean = TRUE, othDepVar = c( colMedians( yMat * 1 ) ), 
   algorithm = GenzBretz(),
   calcVCov = TRUE, returnJacobian = TRUE, dummyVars = NULL )
print( summary( margEffCondObsCovD0 ), digits = 3 )
# now with seemingly only dummy variables
margEffCondObsCovDA <- margEff( estResultBFGS, cond = TRUE,
   atMean = TRUE, othDepVar = c( colMedians( yMat * 1 ) ), 
   algorithm = GenzBretz(),
   calcVCov = TRUE, returnJacobian = TRUE, dummyVars = c( "x1", "x2" ) )
print( summary( margEffCondObsCovDA ), digits = c( 3, 3, 2, 2, 2 ) )

# conditional marginal effects
# (assuming that all other dependent variables are one)
margEffCondOne <- margEff( estResultBFGS, cond = TRUE, othDepVar = 1,
   algorithm = GenzBretz() )
round( margEffCondOne, 3 )

# conditional marginal effects with covariance matrix at sample mean
# (assuming that all other dependent variables are one)
margEffCondOneCov <- margEff( estResultBFGS, cond = TRUE, othDepVar = 1,
   data = as.data.frame( t( colMeans( xMat ) ) ), calcVCov = TRUE,
   algorithm = GenzBretz() )
round( margEffCondOneCov, 3 )
round( attr( margEffCondOneCov, "vcov" ), 2 )
round( drop( attr( margEffCondOneCov, "vcov" ) ), 2 )
print( summary( margEffCondOneCov ), digits = c( 3, 3, 2, 2, 2 ) )
# now with using argument 'atMean'
margEffCondOneCovA <- margEff( estResultBFGS, cond = TRUE, othDepVar = 1,
   atMean = TRUE, algorithm = GenzBretz(), calcVCov = TRUE )
all.equal( margEffCondOneCovA, margEffCondOneCov )

# marginal effects (with factor as explanatory variable, 2 levels)
# unconditional marginal effects (with Jacobian)
dat <- datFac
margEffFacUnc <- margEff( estResultFac, calcVCov = TRUE, returnJacobian = TRUE )
all.equal( margEffUnc, margEffFacUnc, check.attributes = FALSE )
round( margEffFacUnc, 3 )
# now at mean values
margEffFacUncMean <- try( margEff( estResultFac, calcVCov = TRUE, atMean = TRUE ) )

# conditional marginal effects
# (assuming that all other dependent variables are as observed)
margEffFacCondObs <- margEff( estResultFac, cond = TRUE,
   algorithm = GenzBretz() )
all.equal( margEffCondObs, margEffFacCondObs, check.attributes = FALSE )
round( margEffFacCondObs, 3 )
# now at mean values
margEffFacCondObs <- try( margEff( estResultFac, cond = TRUE,
   algorithm = GenzBretz(), atMean = TRUE ) )


# marginal effects (with factor as explanatory variable, 3 levels)
# unconditional marginal effects (with Jacobian)
margEffFac13Unc <- margEff( estResultFac13, calcVCov = TRUE, returnJacobian = TRUE )
round( margEffFac13Unc, 3 )
# now at mean values
margEffFac13UncMean <- try( margEff( estResultFac13, calcVCov = TRUE, atMean = TRUE ) )

# conditional marginal effects
# (assuming that all other dependent variables are as observed)
margEffFac13CondObs <- margEff( estResultFac13, cond = TRUE,
   algorithm = GenzBretz() )
round( margEffFac13CondObs, 3 )
# now at mean values
margEffFac13CondObs <- try( margEff( estResultFac13, cond = TRUE,
   algorithm = GenzBretz(), atMean = TRUE ) )
