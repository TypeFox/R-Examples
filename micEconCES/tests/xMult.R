# load the micEconCES package
library( "micEconCES" )
options( digits = 3 )

# seed for random number generation
set.seed( 123 )

# number of observations
nObs <- 200

# number of explanatory variables
nExog <- 4

# create data set
cesData <- data.frame( obsNo = 1:nObs )

# names of explanatory variables
xxNames <- paste( "xx", 1:nExog, sep = "." )

# add explanatory variables
for( i in 1:nExog ) {
   cesData[[ xxNames[ i ] ]] <- rchisq( nObs, 10 + i )
}
cesData$time <- c( 0:( nObs - 1 ) )

# coefficients
cesCoef <- c( 1, 1:nExog / sum( 1:nExog ), 0.5, 1.1 )
names( cesCoef ) <- c( "gamma", paste( "delta", 1:nExog, sep = "_" ),
   "rho", "nu" )
cesCoefTc <- c( cesCoef[ 1 ], lambda = 0.015, cesCoef[ -1 ] )

# calculate deterministic endogenous variable
cesData$y <- cesCalc( xNames = xxNames, data = cesData, coef = cesCoef )
print( cesData$y )
# check if removing the names of the coefficients makes a difference
all.equal( cesData$y,
   cesCalc( xNames = xxNames, data = cesData, coef = unname( cesCoef ) ) )
# check if permuting the coefficients makes a difference
all.equal( cesData$y,
   cesCalc( xNames = xxNames, data = cesData, coef = sample( cesCoef, 7 ) ) )

# deterministic dependent variable with technological change
cesData$yTc <- cesCalc( xNames = xxNames, tName = "time", data = cesData, 
   coef = cesCoefTc )
print( cesData$yTc )
all.equal( cesData$yTc, 
   cesData$y * exp( cesCoefTc[ "lambda" ] * c( 0:( nObs - 1 ) ) ) )
# check if removing the names of the coefficients makes a difference
all.equal( cesData$yTc,
   cesCalc( xNames = xxNames, tName = "time", data = cesData, 
      coef = unname( cesCoefTc ) ) )
# check if permuting the coefficients makes a difference
all.equal( cesData$yTc,
   cesCalc( xNames = xxNames, tName = "time", data = cesData, 
      coef = sample( cesCoefTc, 8 ) ) )


# adding noise to the endogenous variable
cesData$y <- cesData$y + rnorm( nObs )
