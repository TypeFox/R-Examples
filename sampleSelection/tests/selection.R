library( "sampleSelection" )
library( "mvtnorm" )
library( "lmtest" )
options( digits = 3 )
N <- 1500
NNA <- 5
vc <- diag(3)
vc[lower.tri(vc)] <- c(0.9, 0.5, 0.6)
vc[upper.tri(vc)] <- vc[lower.tri(vc)]
set.seed(1)
## ------- Tobit-5 example ---------
eps <- rmvnorm( N, rep(0, 3), vc )
t5Dat <- data.frame( xs = runif(N) )
t5Dat$ys <- t5Dat$xs + eps[,1] > 0
t5Dat$xo1 <- runif(N)
t5Dat$yo1 <- t5Dat$xo1 + eps[,2]
t5Dat$xo2 <- runif(N)
t5Dat$yo2 <- t5Dat$xo2 + eps[,3]
## Put some NA-s into the data
t5Dat$ys[sample(N, NNA)] <- NA
t5Dat$xs[sample(N, NNA)] <- NA
t5Dat$xo1[sample(N, NNA)] <- NA
t5Dat$xo2[sample(N, NNA)] <- NA
t5Dat$yo1[sample(N, NNA)] <- NA
t5Dat$yo2[sample(N, NNA)] <- NA

## 2-step estimation
testTobit5TwoStep <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2), 
   method = "2step", data = t5Dat )
print( testTobit5TwoStep )
print( summary( testTobit5TwoStep ) )
print( coef( testTobit5TwoStep ) )
print( coef( testTobit5TwoStep, part = "outcome" ) )
print( coef( summary( testTobit5TwoStep ) ) )
print( coef( summary( testTobit5TwoStep ), part = "outcome" ) )
stdEr( testTobit5TwoStep )
print( vcov( testTobit5TwoStep ) )
print( vcov( testTobit5TwoStep, part = "outcome" ) )
nobs( testTobit5TwoStep )
nObs( testTobit5TwoStep )
print( fitted( testTobit5TwoStep, part = "outcome" ) )
print( fitted( testTobit5TwoStep, part = "selection" ) )
print( residuals( testTobit5TwoStep, part = "outcome" ) )
all.equal( residuals( testTobit5TwoStep ),
   residuals( testTobit5TwoStep, part = "outcome" ) )
print( residuals( testTobit5TwoStep, part = "selection" ) )
all.equal( residuals( testTobit5TwoStep, part = "selection" ),
   residuals( testTobit5TwoStep, part = "selection", type = "deviance" ) )
all.equal( residuals( testTobit5TwoStep$probit ),
   residuals( testTobit5TwoStep, part = "selection" ) )
print( residuals( testTobit5TwoStep, part = "selection", type = "pearson" ) )
all.equal( residuals( testTobit5TwoStep$probit, type = "pearson" ),
   residuals( testTobit5TwoStep, part = "selection", type = "pearson" ) )
print( residuals( testTobit5TwoStep, part = "selection", type = "response" ) )
all.equal( residuals( testTobit5TwoStep$probit, type = "response" ),
   residuals( testTobit5TwoStep, part = "selection", type = "response" ) )
t5Samp <- rownames( t5Dat ) %in% names( residuals( testTobit5TwoStep ) )
all.equal( residuals( testTobit5TwoStep, part = "selection", type = "response" ),
   t5Dat$ys[ t5Samp ] - fitted( testTobit5TwoStep, part = "selection" ) )
round( predict( testTobit5TwoStep, newdata = t5Dat, type = "link" ), 3 )
round( predict( testTobit5TwoStep, type = "link" ), 3 )
all.equal(
   predict( testTobit5TwoStep$probit, newdata = t5Dat[ t5Samp, ], type = "link" ),
   predict( testTobit5TwoStep, newdata = t5Dat[ t5Samp, ], type = "link" ) )
all.equal(
   predict( testTobit5TwoStep, newdata = t5Dat[ , "xs", drop = FALSE ],
      type = "link" ),
   predict( testTobit5TwoStep, newdata = t5Dat, type = "link" ) )
round( predict( testTobit5TwoStep, newdata = t5Dat, type = "response" ), 3 )
round( predict( testTobit5TwoStep, type = "response" ), 3 )
all.equal(
   predict( testTobit5TwoStep$probit, newdata = t5Dat[ t5Samp, ], type = "response" ),
   predict( testTobit5TwoStep, newdata = t5Dat[ t5Samp, ], type = "response" ) )
all.equal(
   predict( testTobit5TwoStep, newdata = t5Dat[ , "xs", drop = FALSE ],
      type = "response" ),
   predict( testTobit5TwoStep, newdata = t5Dat, type = "response" ) )
round( predict( testTobit5TwoStep, newdata = t5Dat, type = "unconditional" ), 3 )
all( is.na( predict( testTobit5TwoStep, type = "unconditional" )[
   cbind( t5Dat$ys, !t5Dat$ys )[ t5Samp, ] ] ) )
all.equal( predict( testTobit5TwoStep, type = "unconditional" )[
   !t5Dat$ys[ t5Samp ], 1 ], 
   predict( testTobit5TwoStep, newdata = t5Dat[ t5Samp & !t5Dat$ys, ],
      type = "unconditional" )[ , 1 ] )
all.equal( predict( testTobit5TwoStep, type = "unconditional" )[
   t5Dat$ys[ t5Samp ], 2 ], 
   predict( testTobit5TwoStep, newdata = t5Dat[ t5Samp & t5Dat$ys, ],
      type = "unconditional" )[ , 2 ] )
all.equal(
   predict( testTobit5TwoStep, newdata = t5Dat[ , c( "xo1", "xo2" ) ],
      type = "unconditional" ), 
   predict( testTobit5TwoStep, newdata = t5Dat, type = "unconditional" ) )
round( predict( testTobit5TwoStep, newdata = t5Dat, type = "conditional" ), 3 )
all( is.na( predict( testTobit5TwoStep, type = "conditional" )[
   cbind( t5Dat$ys, t5Dat$ys, !t5Dat$ys, !t5Dat$ys )[ t5Samp, ] ] ) )
all.equal( predict( testTobit5TwoStep, type = "conditional" )[
   !t5Dat$ys[ t5Samp ], 1 ], 
   predict( testTobit5TwoStep, newdata = t5Dat[ t5Samp & !t5Dat$ys, ],
      type = "conditional" )[ , 1 ] )
all.equal( predict( testTobit5TwoStep, type = "conditional" )[
   t5Dat$ys[ t5Samp ], 4 ], 
   predict( testTobit5TwoStep, newdata = t5Dat[ t5Samp & t5Dat$ys, ],
      type = "conditional" )[ , 4 ] )
all.equal(
   rowSums( predict( testTobit5TwoStep, type = "conditional" )[ , c(1,4)], na.rm = TRUE ),
   fitted( testTobit5TwoStep ), check.attributes = FALSE )
all.equal(
   predict( testTobit5TwoStep, newdata = t5Dat[ , c( "xs", "xo1", "xo2" ) ],
      type = "conditional" ), 
   predict( testTobit5TwoStep, newdata = t5Dat, type = "conditional" ) )
mmoTestTobit5TwoStep <- model.matrix( testTobit5TwoStep, part = "outcome" )
print( mmoTestTobit5TwoStep )
mmsTestTobit5TwoStep <- model.matrix( testTobit5TwoStep, part = "selection" )
print( mmsTestTobit5TwoStep )
mfTestTobit5TwoStep <- model.frame( testTobit5TwoStep )
print( mfTestTobit5TwoStep )
try( logLik( testTobit5TwoStep ) )

# the same dependent variables in both of the outcome equations
testTobit5STwoStep <- selection( ys~xs, list(yo1 ~ xo1, yo2 ~ xo1 ),
   method = "2step", data = t5Dat )
print( testTobit5STwoStep )
print( summary( testTobit5STwoStep ) )
nobs( testTobit5STwoStep )
nObs( testTobit5STwoStep )
fitted( testTobit5STwoStep, part = "outcome" )
residuals( testTobit5STwoStep, part = "outcome" )
all.equal( residuals( testTobit5STwoStep ),
   residuals( testTobit5STwoStep, part = "outcome" ) )
t5SSamp <- rownames( t5Dat ) %in% names( residuals( testTobit5STwoStep ) )
round( predict( testTobit5STwoStep, newdata = t5Dat, type = "unconditional" ), 3 )
all( is.na( predict( testTobit5STwoStep, type = "unconditional" )[
   cbind( t5Dat$ys, !t5Dat$ys )[ t5SSamp, ] ] ) )
all.equal( predict( testTobit5STwoStep, type = "unconditional" )[
   !t5Dat$ys[ t5SSamp ], 1 ], 
   predict( testTobit5STwoStep, newdata = t5Dat[ t5SSamp & !t5Dat$ys, ],
      type = "unconditional" )[ , 1 ] )
all.equal( predict( testTobit5STwoStep, type = "unconditional" )[
   t5Dat$ys[ t5SSamp ], 2 ], 
   predict( testTobit5STwoStep, newdata = t5Dat[ t5SSamp & t5Dat$ys, ],
      type = "unconditional" )[ , 2 ] )
all.equal(
   rowSums( predict( testTobit5STwoStep, type = "conditional" )[ , c(1,4)], na.rm = TRUE ),
   fitted( testTobit5STwoStep ), check.attributes = FALSE )
all.equal(
   predict( testTobit5STwoStep, newdata = t5Dat[ , c( "xo1", "xo2" ) ],
      type = "unconditional" ),
   predict( testTobit5STwoStep, newdata = t5Dat, type = "unconditional" ) )
round( predict( testTobit5STwoStep, newdata = t5Dat, type = "conditional" ), 3 )
all( is.na( predict( testTobit5STwoStep, type = "conditional" )[
   cbind( t5Dat$ys, t5Dat$ys, !t5Dat$ys, !t5Dat$ys )[ t5SSamp, ] ] ) )
all.equal( predict( testTobit5STwoStep, type = "conditional" )[
   !t5Dat$ys[ t5SSamp ], 1 ], 
   predict( testTobit5STwoStep, newdata = t5Dat[ t5SSamp & !t5Dat$ys, ],
      type = "conditional" )[ , 1 ] )
all.equal( predict( testTobit5STwoStep, type = "conditional" )[
   t5Dat$ys[ t5SSamp ], 4 ], 
   predict( testTobit5STwoStep, newdata = t5Dat[ t5SSamp & t5Dat$ys, ],
      type = "conditional" )[ , 4 ] )
all.equal(
   predict( testTobit5STwoStep, newdata = t5Dat[ , c( "xs", "xo1", "xo2" ) ],
      type = "conditional" ),
   predict( testTobit5STwoStep, newdata = t5Dat, type = "conditional" ) )
mmsTestTobit5STwoStep <- model.matrix( testTobit5STwoStep, part = "selection" )
print( mmsTestTobit5STwoStep )
mmoTestTobit5STwoStep <- model.matrix( testTobit5STwoStep, part = "outcome" )
print( mmoTestTobit5STwoStep )
mfTestTobit5STwoStep <- model.frame( testTobit5STwoStep )
print( mfTestTobit5STwoStep )

## Tobit-5: ML estimatiom
testTobit5Ml <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2), method = "ml",
   data = t5Dat )
print( testTobit5Ml )
print( summary( testTobit5Ml ) )
print( coef( testTobit5Ml ) )
print( coef( testTobit5Ml, part = "outcome" ) )
print( coef( summary( testTobit5Ml ) ) )
print( coef( summary( testTobit5Ml ), part = "outcome" ) )
stdEr( testTobit5Ml )
print( vcov( testTobit5Ml ) )
print( vcov( testTobit5Ml, part = "outcome" ) )
nobs( testTobit5Ml )
nObs( testTobit5Ml )
print( fitted( testTobit5Ml, part = "outcome" ) )
print( fitted( testTobit5Ml, part = "selection" ) )
print( residuals( testTobit5Ml, part = "outcome" ) )
all.equal( residuals( testTobit5Ml ),
   residuals( testTobit5Ml, part = "outcome" ) )
print( residuals( testTobit5Ml, part = "selection" ) )
all.equal( residuals( testTobit5Ml, part = "selection" ),
   residuals( testTobit5Ml, part = "selection", type = "deviance" ) )
!isTRUE( all.equal( residuals( testTobit5TwoStep, part = "selection" ),
   residuals( testTobit5Ml, part = "selection" ) ) )
print( residuals( testTobit5Ml, part = "selection", type = "pearson" ) )
!isTRUE( all.equal(
   residuals( testTobit5TwoStep, part = "selection", type = "pearson" ),
   residuals( testTobit5Ml, part = "selection", type = "pearson" ) ) )
print( residuals( testTobit5Ml, part = "selection", type = "response" ) )
!isTRUE( all.equal(
   residuals( testTobit5TwoStep, part = "selection", type = "response" ),
   residuals( testTobit5Ml, part = "selection", type = "response" ) ) )
all.equal( residuals( testTobit5Ml, part = "selection", type = "response" ),
   t5Dat$ys[ t5Samp ] - fitted( testTobit5Ml, part = "selection" ) )
all.equal( residuals( testTobit5TwoStep, part = "selection" ),
   residuals( testTobit5TwoStep, part = "selection", type = "deviance" ) )
round( predict( testTobit5Ml, newdata = t5Dat, type = "link" ), 3 )
all.equal( predict( testTobit5Ml, type = "link" ),
   predict( testTobit5Ml, newdata = t5Dat[ t5Samp, ], type = "link" ) )
all.equal(
   predict( testTobit5Ml, newdata = t5Dat[ t5Samp, ], type = "link" ),
   qnorm( fitted( testTobit5Ml, part = "selection" ) ) )
all.equal(
   predict( testTobit5Ml, newdata = t5Dat[ , "xs", drop = FALSE ],
      type = "link" ),
   predict( testTobit5Ml, newdata = t5Dat, type = "link" ) )
round( predict( testTobit5Ml, newdata = t5Dat, type = "response" ), 3 )
all.equal( predict( testTobit5Ml, type = "response" ),
   predict( testTobit5Ml, newdata = t5Dat[ t5Samp, ], type = "response" ) )
all.equal(
   predict( testTobit5Ml, newdata = t5Dat[ t5Samp, ], type = "response" ),
   fitted( testTobit5Ml, part = "selection" ) )
all.equal(
   predict( testTobit5Ml, newdata = t5Dat[ , "xs", drop = FALSE ],
      type = "response" ),
   predict( testTobit5Ml, newdata = t5Dat, type = "response" ) )
round( predict( testTobit5Ml, newdata = t5Dat, type = "unconditional" ), 3 )
all( is.na( predict( testTobit5Ml, type = "unconditional" )[
   cbind( t5Dat$ys, !t5Dat$ys )[ t5Samp, ] ] ) )
all.equal( predict( testTobit5Ml, type = "unconditional" )[
   !t5Dat$ys[ t5Samp ], 1 ], 
   predict( testTobit5Ml, newdata = t5Dat[ t5Samp & !t5Dat$ys, ],
      type = "unconditional" )[ , 1 ] )
all.equal( predict( testTobit5Ml, type = "unconditional" )[
   t5Dat$ys[ t5Samp ], 2 ], 
   predict( testTobit5Ml, newdata = t5Dat[ t5Samp & t5Dat$ys, ],
      type = "unconditional" )[ , 2 ] )
all.equal(
   rowSums( predict( testTobit5Ml, type = "unconditional" ), na.rm = TRUE ),
   fitted( testTobit5Ml ), check.attributes = FALSE )
all.equal(
   predict( testTobit5Ml, newdata = t5Dat[ , c( "xo1", "xo2" ) ],
      type = "unconditional" ),
   predict( testTobit5Ml, newdata = t5Dat, type = "unconditional" ) )
round( predict( testTobit5Ml, newdata = t5Dat, type = "conditional" ), 3 )
all( is.na( predict( testTobit5Ml, type = "conditional" )[
   cbind( t5Dat$ys, t5Dat$ys, !t5Dat$ys, !t5Dat$ys )[ t5Samp, ] ] ) )
all.equal( predict( testTobit5Ml, type = "conditional" )[
   !t5Dat$ys[ t5Samp ], 1 ], 
   predict( testTobit5Ml, newdata = t5Dat[ t5Samp & !t5Dat$ys, ],
      type = "conditional" )[ , 1 ] )
all.equal( predict( testTobit5Ml, type = "conditional" )[
   t5Dat$ys[ t5Samp ], 4 ], 
   predict( testTobit5Ml, newdata = t5Dat[ t5Samp & t5Dat$ys, ],
      type = "conditional" )[ , 4 ] )
all.equal(
   predict( testTobit5Ml, newdata = t5Dat[ , c( "xs", "xo1", "xo2" ) ],
      type = "conditional" ),
   predict( testTobit5Ml, newdata = t5Dat, type = "conditional" ) )
mmsTestTobit5Ml <- model.matrix( testTobit5Ml, part = "selection" )
print( mmsTestTobit5Ml )
mmoTestTobit5Ml <- model.matrix( testTobit5Ml, part = "outcome" )
print( mmoTestTobit5Ml )
mfTestTobit5Ml <- model.frame( testTobit5Ml )
print( mfTestTobit5Ml )
logLik( testTobit5Ml )

# LR tests
testTobit5Ml00 <- selection(ys~xs, list(yo1 ~ 1, yo2 ~ 1), method = "ml",
   data = t5Dat[ t5Samp, ] )
lrtest( testTobit5Ml00, testTobit5Ml )
testTobit5Ml10 <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ 1), method = "ml",
   data = t5Dat[ t5Samp, ] )
lrtest( testTobit5Ml10, testTobit5Ml )
testTobit5Ml01 <- selection(ys~xs, list(yo1 ~ 1, yo2 ~ xo2), method = "ml",
   data = t5Dat[ t5Samp, ] )
lrtest( testTobit5Ml01, testTobit5Ml )

# ML with model.matrices returned
testTobit5MlMm <- selection( ys ~ xs, list( yo1 ~ xo1, yo2 ~ xo2 ),
   method = "ml", xs = TRUE, xo = TRUE, data = t5Dat )
mmsTestTobit5MlMm <- model.matrix( testTobit5MlMm, part = "selection" )
attr( mmsTestTobit5MlMm, "assign" ) <- attr( mmsTestTobit5Ml, "assign" )
all.equal( mmsTestTobit5Ml, mmsTestTobit5MlMm )
mmoTestTobit5MlMm <- model.matrix( testTobit5MlMm, part = "outcome" )
attr( mmoTestTobit5MlMm[[ 1 ]], "assign" ) <- attr( mmoTestTobit5Ml[[ 1 ]], "assign" )
attr( mmoTestTobit5MlMm[[ 2 ]], "assign" ) <- attr( mmoTestTobit5Ml[[ 2 ]], "assign" )
all.equal( mmoTestTobit5Ml, mmoTestTobit5MlMm )
# ML with model.frames returned
testTobit5MlMf <- selection( ys~xs, list( yo1 ~ xo1, yo2 ~ xo2 ),
   method = "ml", mfs = TRUE, mfo = TRUE, data = t5Dat )
mfTestTobit5MlMf <- model.frame( testTobit5MlMf )
all.equal( mfTestTobit5Ml, mfTestTobit5MlMf )

# return just the model.frame
all.equal( selection( ys~xs, list( yo1 ~ xo1, yo2 ~ xo2 ),
   method = "model.frame", data = t5Dat ), mfTestTobit5MlMf )

# the same dependent variables in both of the outcome equations
testTobit5SMl <- selection( ys~xs, list(yo1 ~ xo1, yo2 ~ xo1 ), method = "ml",
   data = t5Dat )
print( testTobit5SMl )
print( summary( testTobit5SMl ) )
nobs( testTobit5SMl )
nObs( testTobit5SMl )
fitted( testTobit5SMl, part = "outcome" )
residuals( testTobit5SMl, part = "outcome" )
all.equal( residuals( testTobit5SMl ),
   residuals( testTobit5SMl, part = "outcome" ) )
round( predict( testTobit5SMl, newdata = t5Dat, type = "unconditional" ), 3 )
all( is.na( predict( testTobit5SMl, type = "unconditional" )[
   cbind( t5Dat$ys, !t5Dat$ys )[ t5SSamp, ] ] ) )
all.equal( predict( testTobit5SMl, type = "unconditional" )[
   !t5Dat$ys[ t5SSamp ], 1 ], 
   predict( testTobit5SMl, newdata = t5Dat[ t5SSamp & !t5Dat$ys, ],
      type = "unconditional" )[ , 1 ] )
all.equal( predict( testTobit5SMl, type = "unconditional" )[
   t5Dat$ys[ t5SSamp ], 2 ], 
   predict( testTobit5SMl, newdata = t5Dat[ t5SSamp & t5Dat$ys, ],
      type = "unconditional" )[ , 2 ] )
all.equal(
   rowSums( predict( testTobit5SMl, type = "unconditional" ), na.rm = TRUE ),
   fitted( testTobit5SMl ), check.attributes = FALSE )
all.equal(
   predict( testTobit5SMl, newdata = t5Dat[ , c( "xo1", "xo2" ) ],
      type = "unconditional" ),
   predict( testTobit5SMl, newdata = t5Dat, type = "unconditional" ) )
round( predict( testTobit5SMl, newdata = t5Dat, type = "conditional" ), 3 )
all( is.na( predict( testTobit5SMl, type = "conditional" )[
   cbind( t5Dat$ys, t5Dat$ys, !t5Dat$ys, !t5Dat$ys )[ t5SSamp, ] ] ) )
all.equal( predict( testTobit5SMl, type = "conditional" )[
   !t5Dat$ys[ t5SSamp ], 1 ], 
   predict( testTobit5SMl, newdata = t5Dat[ t5SSamp & !t5Dat$ys, ],
      type = "conditional" )[ , 1 ] )
all.equal( predict( testTobit5SMl, type = "conditional" )[
   t5Dat$ys[ t5SSamp ], 4 ], 
   predict( testTobit5SMl, newdata = t5Dat[ t5SSamp & t5Dat$ys, ],
      type = "conditional" )[ , 4 ] )
all.equal(
   predict( testTobit5SMl, newdata = t5Dat[ , c( "xs", "xo1", "xo2" ) ],
      type = "conditional" ),
   predict( testTobit5SMl, newdata = t5Dat, type = "conditional" ) )
mmsTestTobit5SMl <- model.matrix( testTobit5SMl, part = "selection" )
print( mmsTestTobit5SMl )
mmoTestTobit5SMl <- model.matrix( testTobit5SMl, part = "outcome" )
print( mmoTestTobit5SMl )
mfTestTobit5SMl <- model.frame( testTobit5SMl )
print( mfTestTobit5SMl )
logLik( testTobit5SMl )

# factors as dependent variable (from Achim Zeileis)
testTobit5FacTwoStep <- selection( factor( ys ) ~ xs,
   list( yo1 ~ xo1, yo2 ~ xo2 ), method = "2step", data = t5Dat )
all.equal( testTobit5FacTwoStep[ -c( 8, 9 ) ], testTobit5TwoStep[ -c( 8, 9 ) ] )
testTobit5YesTwoStep <- selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs,
   list( yo1 ~ xo1, yo2 ~ xo2 ), method = "2step", data = t5Dat )
all.equal( testTobit5YesTwoStep[ -c( 8, 9, 17 ) ],
   testTobit5TwoStep[ -c( 8, 9, 17 ) ] )
all.equal( testTobit5YesTwoStep$param[ -c( 12 ) ],
   testTobit5TwoStep$param[ -c( 12 ) ] )

testTobit5FacMl <- selection( factor( ys ) ~ xs,
   list( yo1 ~ xo1, yo2 ~ xo2 ), data = t5Dat )
all.equal( testTobit5FacMl[ -c( 15, 18, 19 ) ],
   testTobit5Ml[ -c( 15, 18, 19 ) ] )
testTobit5YesMl <- selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs,
   list( yo1 ~ xo1, yo2 ~ xo2 ), data = t5Dat )
all.equal( testTobit5YesMl[ -c( 15, 17, 18, 19 ) ],
   testTobit5Ml[ -c( 15, 17, 18, 19 ) ] )
all.equal( testTobit5YesMl$param[ -c( 10 ) ], testTobit5Ml$param[ -c( 10 ) ] )

# with pre-defined list of outcome equations (works since revision 1420)
oList <- list( yo1 ~ xo1, yo2 ~ xo2 )
testTobit5LiTwoStep <- selection( ys ~ xs, oList, method = "2step",
   data = t5Dat )
all.equal( testTobit5LiTwoStep[ -8 ],  testTobit5TwoStep[ -8 ] )
testTobit5LiFacTwoStep <- selection( factor( ys ) ~ xs, oList, method = "2step",
   data = t5Dat )
all.equal( testTobit5LiFacTwoStep[ -8 ],  testTobit5FacTwoStep[ -8 ] )
testTobit5LiYesTwoStep <- selection( factor( ys, labels = c( "no", "yes" ) )
   ~ xs, oList, method = "2step", data = t5Dat )
all.equal( testTobit5LiYesTwoStep[ -8 ],  testTobit5YesTwoStep[ -8 ] )

testTobit5LiMl <- selection( ys ~ xs, oList, data = t5Dat )
all.equal( testTobit5LiMl[ -18 ],  testTobit5Ml[ -18 ] )
testTobit5LiFacMl <- selection( factor( ys ) ~ xs, oList, data = t5Dat )
all.equal( testTobit5LiFacMl[ -18 ],  testTobit5FacMl[ -18 ] )
testTobit5LiYesMl <- selection( factor( ys, labels = c( "no", "yes" ) )
   ~ xs, oList, data = t5Dat )
all.equal( testTobit5LiYesMl[ -18 ],  testTobit5YesMl[ -18 ] )

# return just the model.frame
selection( factor( ys ) ~ xs, list( yo1 ~ xo1, yo2 ~ xo2 ),
   method = "model.frame", data = t5Dat )
selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs,
   list( yo1 ~ xo1, yo2 ~ xo2 ), method="model.frame", data = t5Dat )

# the case without intercepts 
cat("Now run tobit5 without intercepts\n")
print(coef(selection( ys ~ xs - 1, list( yo1 ~ xo1 - 1, yo2 ~ xo2 - 1),
   data = t5Dat ) ) )
# return just the model.frame
selection( ys ~ xs - 1, list( yo1 ~ xo1 - 1, yo2 ~ xo2 - 1 ),
   method = "model.frame", data = t5Dat )

## estimations withs weights that do not work
testTobit5TwoStepWe <- selection( ys ~ xs, list( yo1 ~ xo1, yo2 ~ xo2), 
   method = "2step", weights = rep( 0.5, N ), data = t5Dat )
all.equal( testTobit5TwoStepWe[-8], testTobit5TwoStep[-8] )

testTobit5MlWe <- selection( ys ~ xs, list( yo1 ~ xo1, yo2 ~ xo2), 
   method = "ml", weights = rep( 0.5, N ), data = t5Dat )
all.equal( testTobit5MlWe[-18], testTobit5Ml[-18] )

## data directly in the workspace
ys <- t5Dat$ys
xs <- t5Dat$xs
yo1 <- t5Dat$yo1
xo1 <- t5Dat$xo1
yo2 <- t5Dat$yo2
xo2 <- t5Dat$xo2

testTobit5WsTwoStep <- selection( ys ~ xs, list( yo1 ~ xo1, yo2 ~ xo2 ), 
   method = "2step" )
all.equal( testTobit5WsTwoStep[-8], testTobit5TwoStep[-8] )
all.equal( model.matrix( testTobit5WsTwoStep ),
   mmoTestTobit5TwoStep )
all.equal( model.matrix( testTobit5WsTwoStep, part = "selection" ),
   mmsTestTobit5TwoStep )
all.equal( model.frame( testTobit5WsTwoStep ),
   mfTestTobit5TwoStep )

testTobit5WsMl <- selection( ys ~ xs, list( yo1 ~ xo1, yo2 ~ xo2 ) )
all.equal( testTobit5Ml[-17], testTobit5Ml[-17] )
all.equal( model.matrix( testTobit5WsMl ),
   mmoTestTobit5Ml )
all.equal( model.matrix( testTobit5WsMl, part = "selection" ),
   mmsTestTobit5Ml )
all.equal( model.frame( testTobit5WsMl ),
   mfTestTobit5Ml )

rm( ys, xs, yo1, xo1, yo2, xo2 )


## ------- Tobit-2 exmple -----------
vc <- diag(2)
vc[2,1] <- vc[1,2] <- -0.7
eps <- rmvnorm( N, rep(0, 2), vc )
t2Dat <- data.frame( xs = runif(N) )
t2Dat$ys <- t2Dat$xs + eps[,1] > 0
t2Dat$xo <- runif(N)
t2Dat$yo <- ( t2Dat$xo + eps[,2])*(t2Dat$ys > 0)
t2Dat$xs[sample(N, NNA)] <- NA
t2Dat$ys[sample(N, NNA)] <- NA
t2Dat$xo[sample(N, NNA)] <- NA
t2Dat$yo[sample(N, NNA)] <- NA
testTobit2TwoStep <- selection( ys ~ xs, yo ~ xo, method = "2step",
   data = t2Dat )
print( testTobit2TwoStep )
print( summary( testTobit2TwoStep ) )
print( coef( testTobit2TwoStep ) )
print( coef( testTobit2TwoStep, part = "outcome" ) )
print( coef( summary( testTobit2TwoStep ) ) )
print( coef( summary( testTobit2TwoStep ), part = "outcome" ) )
stdEr( testTobit2TwoStep )
print( vcov( testTobit2TwoStep ) )
print( vcov( testTobit2TwoStep, part = "outcome" ) )
print( testTobit2TwoStep$invMillsRatio )
nobs( testTobit2TwoStep )
nObs( testTobit2TwoStep )
print( fitted( testTobit2TwoStep, part = "outcome" ) )
all.equal( residuals( testTobit2TwoStep ),
   residuals( testTobit2TwoStep, part = "outcome" ) )
print( fitted( testTobit2TwoStep, part = "selection" ) )
print( residuals( testTobit2TwoStep, part = "outcome" ) )
print( residuals( testTobit2TwoStep, part = "selection" ) )
all.equal( residuals( testTobit2TwoStep, part = "selection" ),
   residuals( testTobit2TwoStep, part = "selection", type = "deviance" ) )
all.equal( residuals( testTobit2TwoStep$probit ),
   residuals( testTobit2TwoStep, part = "selection" ) )
print( residuals( testTobit2TwoStep, part = "selection", type = "pearson" ) )
all.equal( residuals( testTobit2TwoStep$probit, type = "pearson" ),
   residuals( testTobit2TwoStep, part = "selection", type = "pearson" ) )
print( residuals( testTobit2TwoStep, part = "selection", type = "response" ) )
all.equal( residuals( testTobit2TwoStep$probit, type = "response" ),
   residuals( testTobit2TwoStep, part = "selection", type = "response" ) )
t2Samp <- rownames( t2Dat ) %in% names( residuals( testTobit2TwoStep ) )
all.equal( residuals( testTobit2TwoStep, part = "selection", type = "response" ),
   t2Dat$ys[ t2Samp ] - fitted( testTobit2TwoStep, part = "selection" ) )
round( predict( testTobit2TwoStep, newdata = t2Dat, type = "link" ), 3 )
all.equal( predict( testTobit2TwoStep, type = "link" ), 
   predict( testTobit2TwoStep, newdata = t2Dat[ t2Samp, ], type = "link" ) )
all.equal( predict( testTobit2TwoStep$probit, type = "link" ), 
   predict( testTobit2TwoStep, newdata = t2Dat[ t2Samp, ], type = "link" ) )
all.equal(
   predict( testTobit2TwoStep, newdata = t2Dat[ t2Samp, ], type = "link" ),
   qnorm( fitted( testTobit2TwoStep, part = "selection" ) ) )
all.equal(
   predict( testTobit2TwoStep, newdata = t2Dat[ , "xs", drop = FALSE ],
      type = "link" ),
   predict( testTobit2TwoStep, newdata = t2Dat, type = "link" ) )
round( predict( testTobit2TwoStep, newdata = t2Dat, type = "response" ), 3 )
all.equal( predict( testTobit2TwoStep, type = "response" ), 
   predict( testTobit2TwoStep, newdata = t2Dat[ t2Samp, ], type = "response" ) )
all.equal( predict( testTobit2TwoStep$probit, type = "response" ), 
   predict( testTobit2TwoStep, newdata = t2Dat[ t2Samp, ], type = "response" ) )
all.equal(
   predict( testTobit2TwoStep, newdata = t2Dat[ t2Samp, ], type = "response" ),
   fitted( testTobit2TwoStep, part = "selection" ) )
all.equal(
   predict( testTobit2TwoStep, newdata = t2Dat[ , "xs", drop = FALSE ],
      type = "response" ),
   predict( testTobit2TwoStep, newdata = t2Dat, type = "response" ) )
round( predict( testTobit2TwoStep, newdata = t2Dat, type = "unconditional" ), 3 )
all( is.na( predict( testTobit2TwoStep, type = "unconditional" )[
   !t2Dat$ys[ t2Samp ] ] ) )
all.equal( predict( testTobit2TwoStep, type = "unconditional" )[
      t2Dat$ys[ t2Samp ] ], 
   predict( testTobit2TwoStep, newdata = t2Dat[ t2Samp & t2Dat$ys, ],
      type = "unconditional" ) )
all.equal(
   predict( testTobit2TwoStep, newdata = t2Dat[ , "xo", drop = FALSE ],
      type = "unconditional" ),
   predict( testTobit2TwoStep, newdata = t2Dat, type = "unconditional" ) )
round( predict( testTobit2TwoStep, newdata = t2Dat, type = "conditional" ), 3 )
all( is.na( predict( testTobit2TwoStep, type = "conditional" )[
   !t2Dat$ys[ t2Samp ], ] ) )
all.equal( predict( testTobit2TwoStep, type = "conditional" )[
   t2Dat$ys[ t2Samp ], ], 
   predict( testTobit2TwoStep, newdata = t2Dat[ t2Samp & t2Dat$ys, ],
      type = "conditional" ) )
t2oSamp <- !is.na( t2Dat$yo ) & !is.na( t2Dat$xo ) 
all.equal(
   predict( testTobit2TwoStep, newdata = t2Dat[ t2Samp & t2Dat$ys, ],
      type = "conditional" )[ , 2 ],
   fitted( testTobit2TwoStep, part = "outcome" )[ t2Dat$ys[ t2Samp ] ] )
all.equal(
   predict( testTobit2TwoStep, newdata = t2Dat[ , c( "xs", "xo" ) ],
      type = "conditional" ),
   predict( testTobit2TwoStep, newdata = t2Dat, type = "conditional" ) )
mmoTestTobit2TwoStep <- model.matrix( testTobit2TwoStep, part = "outcome" )
print( mmoTestTobit2TwoStep )
mmsTestTobit2TwoStep <- model.matrix( testTobit2TwoStep, part = "selection" )
print( mmsTestTobit2TwoStep )
mfTestTobit2TwoStep <- model.frame( testTobit2TwoStep )
print( mfTestTobit2TwoStep )
try( logLik( testTobit2TwoStep ) )

testTobit2Ml <- selection( ys ~ xs, yo ~ xo, method = "ml", data = t2Dat )
print( testTobit2Ml )
print( summary( testTobit2Ml ) )
print( coef( testTobit2Ml ) )
print( coef( testTobit2Ml, part = "outcome" ) )
print( coef( summary( testTobit2Ml ) ) )
print( coef( summary( testTobit2Ml ), part = "outcome" ) )
stdEr( testTobit2Ml )
print( vcov( testTobit2Ml ) )
print( vcov( testTobit2Ml, part = "outcome" ) )
nobs( testTobit2Ml )
nObs( testTobit2Ml )
print( fitted( testTobit2Ml, part = "outcome" ) )
print( fitted( testTobit2Ml, part = "selection" ) )
print( residuals( testTobit2Ml, part = "outcome" ) )
all.equal( residuals( testTobit2Ml ),
   residuals( testTobit2Ml, part = "outcome" ) )
print( residuals( testTobit2Ml, part = "selection" ) )
all.equal( residuals( testTobit2Ml, part = "selection" ),
   residuals( testTobit2Ml, part = "selection", type = "deviance" ) )
!isTRUE( all.equal( residuals( testTobit2TwoStep, part = "selection" ),
   residuals( testTobit2Ml, part = "selection" ) ) )
print( residuals( testTobit2Ml, part = "selection", type = "pearson" ) )
!isTRUE( all.equal(
   residuals( testTobit2TwoStep, part = "selection", type = "pearson" ),
   residuals( testTobit2Ml, part = "selection", type = "pearson" ) ) )
print( residuals( testTobit2Ml, part = "selection", type = "response" ) )
!isTRUE( all.equal(
   residuals( testTobit2TwoStep, part = "selection", type = "response" ),
   residuals( testTobit2Ml, part = "selection", type = "response" ) ) )
all.equal( residuals( testTobit2Ml, part = "selection", type = "response" ),
   t2Dat$ys[ t2Samp ] - fitted( testTobit2Ml, part = "selection" ) )
round( predict( testTobit2Ml, newdata = t2Dat, type = "link" ), 3 )
all.equal( predict( testTobit2Ml, type = "link" ), 
   predict( testTobit2Ml, newdata = t2Dat[ t2Samp, ], type = "link" ) )
all.equal(
   predict( testTobit2Ml, newdata = t2Dat[ t2Samp, ], type = "link" ),
   qnorm( fitted( testTobit2Ml, part = "selection" ) ) )
all.equal(
   predict( testTobit2Ml, newdata = t2Dat[ , "xs", drop = FALSE ],
      type = "link" ),
   predict( testTobit2Ml, newdata = t2Dat, type = "link" ) )
round( predict( testTobit2Ml, newdata = t2Dat, type = "response" ), 3 )
all.equal( predict( testTobit2Ml, type = "response" ), 
   predict( testTobit2Ml, newdata = t2Dat[ t2Samp, ], type = "response" ) )
all.equal(
   predict( testTobit2Ml, newdata = t2Dat[ t2Samp, ], type = "response" ),
   fitted( testTobit2Ml, part = "selection" ) )
all.equal(
   predict( testTobit2Ml, newdata = t2Dat[ , "xs", drop = FALSE ],
      type = "response" ),
   predict( testTobit2Ml, newdata = t2Dat, type = "response" ) )
round( predict( testTobit2Ml, newdata = t2Dat, type = "unconditional" ), 3 )
all( is.na( predict( testTobit2Ml, type = "unconditional" )[
   !t2Dat$ys[ t2Samp ] ] ) )
all.equal( predict( testTobit2Ml, type = "unconditional" )[
   t2Dat$ys[ t2Samp ] ], 
   predict( testTobit2Ml, newdata = t2Dat[ t2Samp & t2Dat$ys, ],
      type = "unconditional" ) )
all.equal( predict( testTobit2Ml,
   newdata = t2Dat[ t2Samp & t2Dat$ys, ], type = "unconditional" ),
   fitted( testTobit2Ml, part = "outcome" )[ t2Dat$ys[ t2Samp ] ] )
all.equal(
   predict( testTobit2Ml, newdata = t2Dat[ , "xo", drop = FALSE ],
      type = "unconditional" ),
   predict( testTobit2Ml, newdata = t2Dat, type = "unconditional" ) )
round( predict( testTobit2Ml, newdata = t2Dat, type = "conditional" ), 3 )
all( is.na( predict( testTobit2Ml, type = "conditional" )[
   !t2Dat$ys[ t2Samp ] ] ) )
all.equal( predict( testTobit2Ml, type = "conditional" )[
   t2Dat$ys[ t2Samp ], ], 
   predict( testTobit2Ml, newdata = t2Dat[ t2Samp & t2Dat$ys, ],
      type = "conditional" ) )
all.equal(
   predict( testTobit2Ml, newdata = t2Dat[ , c( "xs", "xo" ) ],
      type = "conditional" ),
   predict( testTobit2Ml, newdata = t2Dat, type = "conditional" ) )
mmsTestTobit2Ml <- model.matrix( testTobit2Ml, part = "selection" )
print( mmsTestTobit2Ml )
mmoTestTobit2Ml <- model.matrix( testTobit2Ml, part = "outcome" )
print( mmoTestTobit2Ml )
mfTestTobit2Ml <- model.frame( testTobit2Ml )
print( mfTestTobit2Ml )
logLik( testTobit2Ml )

# LR test
testTobit2Ml0 <- selection( ys ~ xs, yo ~ 1, method = "ml",
   data = t2Dat[ t2Samp, ] )
lrtest( testTobit2Ml0, testTobit2Ml )

# ML with model.matrices returned
testTobit2MlMm <- selection( ys ~ xs, yo ~ xo, method = "ml", 
   xs = TRUE, xo = TRUE, data = t2Dat )
mmsTestTobit2MlMm <- model.matrix( testTobit2MlMm, part = "selection" )
attr( mmsTestTobit2MlMm, "assign" ) <- attr( mmsTestTobit2Ml, "assign" )
all.equal( mmsTestTobit2Ml, mmsTestTobit2MlMm )
mmoTestTobit2MlMm <- model.matrix( testTobit2MlMm, part = "outcome" )
attr( mmoTestTobit2MlMm, "assign" ) <- attr( mmoTestTobit2Ml, "assign" )
all.equal( mmoTestTobit2Ml, mmoTestTobit2MlMm )
# ML with model.frames returned
testTobit2MlMf <- selection( ys ~ xs, yo ~ xo, method = "ml",
   mfs = TRUE, mfo = TRUE, data = t2Dat )
mfTestTobit2MlMf <- model.frame( testTobit2MlMf )
all.equal( mfTestTobit2Ml, mfTestTobit2MlMf )
                           # attributes (terms) differ here, I don't exactly know how to improve that

# return just the model.frame
all.equal( selection( ys ~ xs, yo ~ xo, method = "model.frame", data = t2Dat ),
   mfTestTobit2MlMf )


## two-step estimation with equal weights
t2Dat$we <- rep( 0.7, N )
testTobit2TwoStepWe <- selection( ys ~ xs, yo ~ xo, method = "2step",
   weights = t2Dat$we, data = t2Dat )
summary( testTobit2TwoStepWe )
all.equal( coef( testTobit2TwoStepWe ), coef( testTobit2TwoStep ) )
nobs( testTobit2TwoStepWe )
nObs( testTobit2TwoStepWe )
try( logLik( testTobit2TwoStepWe ) )

## ML estimation with equal weights
testTobit2MlWe <- selection( ys ~ xs, yo ~ xo, weights = t2Dat$we,
   data = t2Dat )
summary( testTobit2MlWe )
all.equal( coef( testTobit2MlWe ), coef( testTobit2Ml ), tol = 1e-4 )
nobs( testTobit2MlWe )
nObs( testTobit2MlWe )
logLik( testTobit2MlWe )

# LR test
testTobit2MlWe0 <- selection( ys ~ xs, yo ~ 1, weights = t2Dat$we[ t2Samp ],
   method = "ml", data = t2Dat[ t2Samp, ] )
lrtest( testTobit2MlWe0, testTobit2MlWe )

## two-step estimation with unequal weights
t2Dat$wu <- 2 * runif( N )
testTobit2TwoStepWu <- selection( ys ~ xs, yo ~ xo, method = "2step",
   weights = t2Dat$wu, data = t2Dat )
summary( testTobit2TwoStepWu )
nobs( testTobit2TwoStepWu )
nObs( testTobit2TwoStepWu )
try( logLik( testTobit2TwoStepWu ) )

## ML estimation with unequal weights
testTobit2MlWu <- selection( ys ~ xs, yo ~ xo, weights = t2Dat$wu,
   data = t2Dat )
summary( testTobit2MlWu )
nobs( testTobit2MlWu )
nObs( testTobit2MlWu )
logLik( testTobit2MlWu )

# LR test
testTobit2MlWu0 <- selection( ys ~ xs, yo ~ 1, weights = t2Dat$wu[ t2Samp ],
   method = "ml", data = t2Dat[ t2Samp, ] )
lrtest( testTobit2MlWu0, testTobit2MlWu )

## estimations with weights that do not work
try( selection( ys ~ xs, yo ~ xo, method = "2step", weights = 1:99,
   data = t2Dat ) )

try( selection( ys ~ xs, yo ~ xo, method = "ml", weights = 4:14,
   data = t2Dat ) )


# factors as dependent variable (from Achim Zeileis)
testTobit2FacTwoStep <- selection( factor( ys ) ~ xs, yo ~ xo,
   method = "2step", data = t2Dat )
all.equal( testTobit2FacTwoStep[-1], testTobit2TwoStep[-1] )
all.equal( fitted( testTobit2FacTwoStep, part = "selection" ),
   fitted( testTobit2TwoStep, part = "selection" ) )
all.equal(
   residuals( testTobit2FacTwoStep, part = "selection", type = "deviance" ),
   residuals( testTobit2TwoStep, part = "selection", type = "deviance" ) )
all.equal(
   residuals( testTobit2FacTwoStep, part = "selection", type = "pearson" ),
   residuals( testTobit2TwoStep, part = "selection", type = "pearson" ) )
all.equal(
   residuals( testTobit2FacTwoStep, part = "selection", type = "response" ),
   residuals( testTobit2TwoStep, part = "selection", type = "response" ) )

testTobit2YesTwoStep <- selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs,
   yo ~ xo, method = "2step", data = t2Dat )
all.equal( testTobit2YesTwoStep[-1], testTobit2TwoStep[-1] )
all.equal( fitted( testTobit2YesTwoStep, part = "selection" ),
   fitted( testTobit2TwoStep, part = "selection" ) )
all.equal(
   residuals( testTobit2YesTwoStep, part = "selection", type = "deviance" ),
   residuals( testTobit2TwoStep, part = "selection", type = "deviance" ) )
all.equal(
   residuals( testTobit2YesTwoStep, part = "selection", type = "pearson" ),
   residuals( testTobit2TwoStep, part = "selection", type = "pearson" ) )
all.equal(
   residuals( testTobit2YesTwoStep, part = "selection", type = "response" ),
   residuals( testTobit2TwoStep, part = "selection", type = "response" ) )

testTobit2FacMl <- selection( factor( ys ) ~ xs, yo ~ xo, data = t2Dat )
all.equal( testTobit2FacMl[-c(18,19)], testTobit2Ml[-c(18,19)] )
all.equal( fitted( testTobit2FacMl, part = "selection" ),
   fitted( testTobit2Ml, part = "selection" ) )
all.equal(
   residuals( testTobit2FacMl, part = "selection", type = "deviance" ),
   residuals( testTobit2Ml, part = "selection", type = "deviance" ) )
all.equal(
   residuals( testTobit2FacMl, part = "selection", type = "pearson" ),
   residuals( testTobit2Ml, part = "selection", type = "pearson" ) )
all.equal(
   residuals( testTobit2FacMl, part = "selection", type = "response" ),
   residuals( testTobit2Ml, part = "selection", type = "response" ) )

testTobit2YesMl <- selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs,
   yo ~ xo, data = t2Dat )
all.equal( testTobit2YesMl[-c(17:19)], testTobit2Ml[-c(17:19)] )
all.equal( testTobit2YesMl$param[-9], testTobit2Ml$param[-9] )
all.equal( fitted( testTobit2YesMl, part = "selection" ),
   fitted( testTobit2Ml, part = "selection" ) )
all.equal(
   residuals( testTobit2YesMl, part = "selection", type = "deviance" ),
   residuals( testTobit2Ml, part = "selection", type = "deviance" ) )
all.equal(
   residuals( testTobit2YesMl, part = "selection", type = "pearson" ),
   residuals( testTobit2Ml, part = "selection", type = "pearson" ) )
all.equal(
   residuals( testTobit2YesMl, part = "selection", type = "response" ),
   residuals( testTobit2Ml, part = "selection", type = "response" ) )

# return just the model.frame
selection( factor( ys ) ~ xs, yo ~ xo, method = "model.frame", data = t2Dat )
selection( factor( ys, labels = c( "no", "yes" ) ) ~ xs, yo ~ xo,
   method = "model.frame", data = t2Dat )

# the case without intercepts (by Lucas Salazar)
cat("Now run tobit2 without intercepts\n")
print( coef( selection( ys ~ xs - 1, yo ~ xo - 1, data = t2Dat ) ) )
# return just the model.frame
selection( ys ~ xs - 1, yo ~ xo - 1, method = "model.frame", data = t2Dat )

## Raphael Abiry: does 'start' argument work?
init <- coef(testTobit2Ml)
testTobit2MlStart <- selection( ys ~ xs, yo ~ xo, data = t2Dat, method = "ml",
   start = init )
print( summary( testTobit2MlStart ) )
                           # Note: should be only 1 iteration
all.equal( testTobit2Ml[ -c(2,5,6,9,15,16,18)],
   testTobit2MlStart[ -c(2,5,6,9,15,16,18)], tol = 1e-5 )

## Chris Hane: dummy variable (factor) as explanatory variable
t2Dat$xF <- rbinom( nrow( t2Dat ), 1, 0.5 )
testTobit2XsFacTwoStep <- selection( ys ~ xs + factor(xF), yo ~ xo,
   data = t2Dat, method = "2step" )
print( testTobit2XsFacTwoStep )
fitted( testTobit2XsFacTwoStep, part = "selection" )

testTobit2XsFacMl <- selection( ys ~ xs + factor(xF), yo ~ xo, data = t2Dat )
print( testTobit2XsFacMl )
fitted( testTobit2XsFacMl, part = "selection" )
logLik( testTobit2XsFacMl )
lrtest( testTobit2XsFacMl, testTobit2Ml )

testTobit2XoFacTwoStep <- selection( ys ~ xs, yo ~ xo + factor(xF),
   data = t2Dat, method = "2step" )
print( testTobit2XoFacTwoStep )
fitted( testTobit2XoFacTwoStep, part = "selection" )

testTobit2XoFacMl <- selection( ys ~ xs, yo ~ xo + factor(xF), data = t2Dat )
print( testTobit2XoFacMl )
fitted( testTobit2XoFacMl, part = "selection" )
logLik( testTobit2XoFacMl )
lrtest( testTobit2XoFacMl, testTobit2Ml )

## data directly in the workspace
ys <- t2Dat$ys
xs <- t2Dat$xs
yo <- t2Dat$yo
xo <- t2Dat$xo

testTobit2WsTwoStep <- selection( ys ~ xs, yo ~ xo, method = "2step" )
all.equal( testTobit2WsTwoStep[-1], testTobit2TwoStep[-1] )
all.equal( model.matrix( testTobit2WsTwoStep ),
   mmoTestTobit2TwoStep )
all.equal( model.matrix( testTobit2WsTwoStep, part = "selection" ),
   mmsTestTobit2TwoStep )
all.equal( model.frame( testTobit2WsTwoStep ),
   mfTestTobit2TwoStep )

testTobit2WsMl <- selection( ys ~ xs, yo ~ xo )
all.equal( testTobit2Ml[-17], testTobit2Ml[-17] )
all.equal( model.matrix( testTobit2WsMl ),
   mmoTestTobit2Ml )
all.equal( model.matrix( testTobit2WsMl, part = "selection" ),
   mmsTestTobit2Ml )
all.equal( model.frame( testTobit2WsMl ),
   mfTestTobit2Ml )

rm( ys, xs, yo, xo )
