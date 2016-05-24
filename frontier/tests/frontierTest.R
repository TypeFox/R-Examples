library( "frontier" )
library( "plm" )
options( digits = 5 )

printAll <- function( x ) {
   for( n in names( x ) ) {
      cat( "$", n, "\n", sep = "" )
      if( n %in% c( "olsParam", "gridAdj", "gridParam", "mleParam", "olsStdEr", 
            "mleCov", "resid", "fitted", "olsResid" ) ) {
         print( round( x[[ n ]], 2 ) )
      } else {
         print( x[[ n ]] )
      }
      cat( "\n" )
   }
   cat( "class\n" )
   print( class( x ) )
}

printME <- function( x ) {
   me <- attr( x, "margEff" )
   attr( x, "margEff" ) <- NULL
   print( round( x, 2 ) )
   if( !is.null( me ) ) {
      cat( "margEff\n" )
      print( round( me, 2 ) )
   }
}

## example data included in FRONTIER 4.1 (cross-section data)
data( front41Data )
row.names( front41Data ) <- paste( "F", row.names( front41Data ), sep = "_" )
front41Data$logOutput  <- log( front41Data$output )
front41Data$logCapital <- log( front41Data$capital )
front41Data$logLabour  <- log( front41Data$labour )
front41Data$firmNo     <- c( 1:nrow( front41Data ) )
front41Data$ones       <- 1
   
## cross-section data, error components frontier
sa1 <- sfa( logOutput ~ logCapital + logLabour, data = front41Data )
Sa1 <- sfa( log( output ) ~ log( capital ) + log( labour ), data = front41Data )
all.equal( Sa1[-42], sa1[-42], check.attributes = FALSE, tol = 1e-4 )
a1 <- frontier( data = front41Data, "logOutput",
   c( "logCapital", "logLabour" ) )
all.equal( sa1[-42], a1[-42], tol = 1e-4 )
sa1i <- sfa( logOutput ~ ones + logCapital + logLabour - 1, data = front41Data )
all.equal( sa1i[ -c( 3, 7, 20, 42 ) ], sa1[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( a1 )
coef( a1, which = "start" )
round( coef( a1, which = "ols" ), 2 )
coef( a1, which = "ols", extraPar = TRUE )
round( coef( a1, which = "grid" ), 2 )
round( coef( a1, which = "grid", extraPar = TRUE ), 2 )
round( coef( a1 ), 2 )
round( coef( a1, extraPar = TRUE ), 2 )
round( coef( summary( a1 ), which = "ols" ), 2 )
round( coef( summary( a1 ) ), 2 )
round( coef( summary( a1, extraPar = TRUE ) ), 2 )
round( vcov( a1 ), 2 )
round( vcov( a1, extraPar = TRUE ), 2 )
print( logLik( a1, which = "ols" ), digits = 4 )
print( logLik( a1, which = "grid" ), digits = 4 )
print( logLik( a1 ), digits = 4 )
nobs( a1 )
print( summary( a1 ), digits = 1 )
print( summary( a1, effMinusU = FALSE ), digits = 1 )
print( summary( a1, extraPar = TRUE ), digits = 1 )
all.equal( summary( sa1, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sa1i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( sa1, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sa1i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( a1 )
printME( efficiencies( a1, margEff = TRUE ) )
printME( efficiencies( a1, asInData = TRUE ) )
printME( efficiencies( a1, minusU = FALSE ) )
printME( efficiencies( a1, asInData = TRUE, minusU = FALSE ) )
all.equal( efficiencies( a1 ), efficiencies( sa1i ) )
all.equal( efficiencies( a1, asInData = TRUE, minusU = FALSE ), 
   efficiencies( sa1i, asInData = TRUE, minusU = FALSE ) )
round( fitted( a1 ), 2 )
round( fitted( a1, asInData = TRUE ), 2 )
round( residuals( a1 ), 2 )
round( residuals( a1, asInData = TRUE ), 2 )
all.equal( fitted( a1, asInData = TRUE ) + residuals( a1, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )
printAll( a1 )

## cross-section data, error components frontier, truncNorm
sa2 <- sfa( logOutput ~ logCapital + logLabour, data = front41Data,
   truncNorm = TRUE, printIter = 4 )
a2 <- frontier( data = front41Data, "logOutput",
   c( "logCapital", "logLabour" ), truncNorm = TRUE )
all.equal( sa2[-c(8,42)], a2[-c(8,42)], tol = 1e-4 )
sa2i <- sfa( logOutput ~ ones + logCapital + logLabour - 1, data = front41Data,
   truncNorm = TRUE )
all.equal( sa2i[ -c( 3, 7, 8, 20, 42 ) ], sa2[ -c( 3, 7, 8, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( a2, digits = 1 )
coef( a2, which = "start" )
round( coef( a2, which = "ols" ), 2 )
round( coef( a2, which = "grid" ), 2 )
round( coef( a2, which = "grid", extraPar = TRUE ), 2 )
round( coef( a2 ), 2 )
round( coef( a2, extraPar = TRUE ), 2 )
round( coef( summary( a2 ), which = "ols" ), 2 )
round( coef( summary( a2 ) ), 2 )
round( coef( summary( a2, extraPar = TRUE ) ), 2 )
round( vcov( a2 ), 2 )
round( vcov( a2, extraPar = TRUE ), 2 )
print( logLik( a2, which = "ols" ), digits = 4 )
print( logLik( a2 ), digits = 4 )
nobs( a2 )
print( summary( a2 ), digits = 1 )
print( summary( a2, extraPar = TRUE ), digits = 1 )
all.equal( summary( sa2, extraPar = TRUE )[ -c( 3, 7, 8, 20, 42 ) ], 
   summary( sa2i, extraPar = TRUE )[ -c( 3, 7, 8, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( sa2, effMinusU = FALSE )[ -c( 3, 7, 8, 20, 42 ) ], 
   summary( sa2i, effMinusU = FALSE )[ -c( 3, 7, 8, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( a2 )
round( efficiencies( a2 ), 2 )
round( efficiencies( a2, asInData = TRUE ), 2 )
all.equal( efficiencies( a2, minusU = FALSE ), 
   efficiencies( sa2i, minusU = FALSE ) )
all.equal( efficiencies( a2, asInData = TRUE ), 
   efficiencies( sa2i, asInData = TRUE ) )
round( residuals( a2 ), 2 )
round( residuals( a2, asInData = TRUE ), 2 )
all.equal( fitted( a2, asInData = TRUE ) + residuals( a2, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )
printAll( a2 )

## cross-section data, error components frontier, truncNorm, starting values
sa5 <- sfa( logOutput ~ logCapital + logLabour, data = front41Data,
   truncNorm = TRUE, startVal = c( 0.5, 0.3, 0.5, 0.5, 0.9, -1 ) )
a5 <- frontier( data = front41Data, "logOutput",
   c( "logCapital", "logLabour" ), truncNorm = TRUE,
   startVal = c( 0.5, 0.3, 0.5, 0.5, 0.9, -1 ) )
all.equal( sa5[-42], a5[-42], tol = 1e-4 )
sa5i <- sfa( logOutput ~ ones + logCapital + logLabour - 1, data = front41Data,
   truncNorm = TRUE, startVal = c( 0.5, 0.3, 0.5, 0.5, 0.9, -1 ) )
all.equal( sa5i[ -c( 3, 7, 21, 42 ) ], sa5[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( a5, digits = 1 )
coef( a5, which = "start" )
round( coef( a5, which = "ols" ), 2 )
coef( a5, which = "grid" )
coef( a5, which = "grid", extraPar = TRUE )
round( coef( a5 ), 2 )
round( coef( a5, extraPar = TRUE ), 2 )
round( coef( summary( a5 ), which = "ols" ), 2 )
round( coef( summary( a5 ) ), 2 )
round( coef( summary( a5, extraPar = TRUE ) ), 2 )
round( vcov( a5 ), 2 )
round( vcov( a5, extraPar = TRUE ), 2 )
print( logLik( a5, which = "ols" ), digits = 4 )
print( logLik( a5 ), digits = 4 )
nobs( a5 )
print( summary( a5 ), digits = 1 )
print( summary( a5, extraPar = TRUE ), digits = 1 )
all.equal( summary( sa5, extraPar = TRUE )[ -c( 3, 7, 21, 42 ) ], 
   summary( sa5i, extraPar = TRUE )[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( sa5, effMinusU = FALSE )[ -c( 3, 7, 21, 42 ) ], 
   summary( sa5i, effMinusU = FALSE )[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE )
lrtest( a5 )
round( efficiencies( a5 ), 2 )
round( efficiencies( a5, asInData = TRUE ), 2 )
all.equal( efficiencies( a5 ), efficiencies( sa5i ) )
all.equal( efficiencies( a5, asInData = TRUE, minusU = FALSE ), 
   efficiencies( sa5i, asInData = TRUE, minusU = FALSE ) )
all.equal( fitted( a5, asInData = TRUE ) + residuals( a5, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )
printAll( a5 )

## cross-section data, efficiency effects frontier
saa1 <- sfa( logOutput ~ logCapital + logLabour | firmNo - 1,
   data = front41Data )
Saa1 <- sfa( log( output ) ~ log( capital ) + log( labour ) | firmNo - 1,
   data = front41Data, printIter = 3 )
all.equal( Saa1[-c(8,42)], saa1[-c(8,42)], check.attributes = FALSE, tol = 1e-4 )
aa1 <- frontier( data = front41Data, "logOutput",
   c( "logCapital", "logLabour" ), zNames = "firmNo" )
all.equal( saa1[-42], aa1[-42], tol = 1e-4 )
saa1i <- sfa( logOutput ~ ones + logCapital + logLabour - 1 | firmNo - 1, 
   data = front41Data )
all.equal( saa1i[ -c( 3, 7, 20, 42 ) ], saa1[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( aa1, digits = 1 )
coef( aa1, which = "start" )
round( coef( aa1, which = "ols" ), 2 )
round( coef( aa1, which = "grid" ), 2 )
round( coef( aa1, which = "grid", extraPar = TRUE ), 2 )
round( coef( aa1 ), 2 )
round( coef( aa1, extraPar = TRUE ), 2 )
round( coef( summary( aa1 ), which = "ols" ), 2 )
round( coef( summary( aa1 ) ), 2 )
round( coef( summary( aa1, extraPar = TRUE ) ), 2 )
round( vcov( aa1 ), 2 )
round( vcov( aa1, extraPar = TRUE ), 2 )
nobs( aa1 )
print( summary( aa1 ), digits = 1 )
print( summary( aa1, effMinusU = FALSE ), digits = 1 )
print( summary( aa1, extraPar = TRUE ), digits = 1 )
all.equal( summary( saa1, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( saa1i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( saa1, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( saa1i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( aa1 )
printME( aa1eff <- efficiencies( aa1, margEff = TRUE ) )
printME( aa1effD <- efficiencies( aa1, asInData = TRUE, margEff = TRUE ) )
printME( aa1effF <- efficiencies( aa1, minusU = FALSE, margEff = TRUE ) )
printME( aa1effDF <- efficiencies( aa1, asInData = TRUE, minusU = FALSE, 
   margEff = TRUE ) )
all.equal( efficiencies( saa1, margEff = TRUE ), 
   efficiencies( saa1i, margEff = TRUE ) )
all.equal( efficiencies( saa1, asInData = TRUE, margEff = TRUE, minusU = FALSE ), 
   efficiencies( saa1i, asInData = TRUE, margEff = TRUE, minusU = FALSE ) )
aa1m <- aa1
aa1m$dataTable[ , "firmNo" ] <- aa1m$dataTable[ , "firmNo" ] + 1e-6
all.equal( attr( aa1eff, "margEff" )[ , 1, 1 ], 
   ( efficiencies( aa1m ) - aa1eff )[ , 1 ] / 1e-6, tol = 1e-4 )
all.equal( attr( aa1effD, "margEff" )[ , 1 ], 
   c( efficiencies( aa1m, asInData = TRUE ) - aa1effD ) / 1e-6, tol = 1e-4 )
all.equal( attr( aa1effF, "margEff" )[ , 1, 1 ], 
   ( efficiencies( aa1m, minusU = FALSE ) - aa1effF )[ , 1 ] / 1e-6, tol = 1e-4 )
all.equal( attr( aa1effDF, "margEff" )[ , 1 ],
   c( efficiencies( aa1m, asInData = TRUE, minusU = FALSE ) - aa1effDF ) / 1e-6,
   tol = 1e-4 )
round( residuals( aa1 ), 2 )
round( residuals( aa1, asInData = TRUE ), 2 )
all.equal( fitted( aa1, asInData = TRUE ) + residuals( aa1, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )
printAll( aa1 )

## cross-section data, efficiency effects frontier, zIntercept
saa2 <- sfa( logOutput ~ logCapital + logLabour | firmNo,
   data = front41Data )
aa2 <- frontier( data = front41Data, "logOutput",
   c( "logCapital", "logLabour" ), zNames = "firmNo", zIntercept = TRUE )
all.equal( saa2[-42], aa2[-42], tol = 1e-4 )
saa2i <- sfa( logOutput ~ ones + logCapital + logLabour - 1 | firmNo,
   data = front41Data )
all.equal( saa2i[ -c( 3, 7, 20, 42 ) ], saa2[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( aa2, digits = 2 )
coef( aa2, which = "start" )
round( coef( aa2, which = "ols" ), 2 )
round( coef( aa2, which = "grid" ), 2 )
round( coef( aa2, which = "grid", extraPar = TRUE ), 2 )
round( coef( aa2 ), 2 )
round( coef( aa2, extraPar = TRUE ), 2 )
round( coef( summary( aa2 ), which = "ols" ), 2 )
round( coef( summary( aa2 ) ), 2 )
round( coef( summary( aa2, extraPar = TRUE ) ), 2 )
round( vcov( aa2 ), 2 )
round( vcov( aa2, extraPar = TRUE ), 2 )
nobs( aa2 )
print( summary( aa2 ), digits = 1 )
print( summary( aa2, extraPar = TRUE ), digits = 1 )
all.equal( summary( saa2, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( saa2i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( saa2, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( saa2i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( aa2 )
printME( efficiencies( aa2, margEff = TRUE ) )
printME( efficiencies( aa2, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( saa2, margEff = TRUE, minusU = FALSE ), 
   efficiencies( saa2i, margEff = TRUE, minusU = FALSE ) )
all.equal( efficiencies( saa2, asInData = TRUE, margEff = TRUE ), 
   efficiencies( saa2i, asInData = TRUE, margEff = TRUE ) )
round( fitted( aa2 ), 2 )
round( fitted( aa2, asInData = TRUE ), 2 )
round( residuals( aa2 ), 2 )
round( residuals( aa2, asInData = TRUE ), 2 )
all.equal( fitted( aa2, asInData = TRUE ) + residuals( aa2, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )
printAll( aa2 )

## cross-section data, efficiency effects frontier, zIntercept, starting values
saa5 <- sfa( logOutput ~ logCapital + logLabour | firmNo,
   data = front41Data, startVal = c( 0.5, 0.3, 0.5, -0.4, -0.01 , 0.4, 0.9 ) )
aa5 <- frontier( data = front41Data, "logOutput",
   c( "logCapital", "logLabour" ), zNames = "firmNo", zIntercept = TRUE,
   startVal = c( 0.5, 0.3, 0.5, -0.4, -0.01 , 0.4, 0.9 ) )
all.equal( saa5[-42], aa5[-42], tol = 1e-4 )
saa5i <- sfa( logOutput ~ ones + logCapital + logLabour - 1 | firmNo, 
   data = front41Data, startVal = c( 0.5, 0.3, 0.5, -0.4, -0.01 , 0.4, 0.9 ) )
all.equal( saa5i[ -c( 3, 7, 21, 42 ) ], saa5[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( aa5, digits = 2 )
coef( aa5, which = "start" )
round( coef( aa5, which = "ols" ), 2 )
coef( aa5, which = "grid" )
coef( aa5, which = "grid", extraPar = TRUE )
round( coef( aa5 ), 2 )
round( coef( aa5, extraPar = TRUE ), 2 )
round( coef( summary( aa5 ), which = "ols" ), 2 )
round( coef( summary( aa5 ) ), 2 )
round( coef( summary( aa5, extraPar = TRUE ) ), 2 )
round( vcov( aa5 ), 2 )
round( vcov( aa5, extraPar = TRUE ), 2 )
nobs( aa5 )
print( summary( aa5 ), digits = 1 )
print( summary( aa5, extraPar = TRUE ), digits = 1 )
all.equal( summary( saa5, extraPar = TRUE )[ -c( 3, 7, 21, 42 ) ], 
   summary( saa5i, extraPar = TRUE )[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( saa5, effMinusU = FALSE )[ -c( 3, 7, 21, 42 ) ], 
   summary( saa5i, effMinusU = FALSE )[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE )
lrtest( aa5 )
printME( efficiencies( aa5, margEff = TRUE ) )
printME( efficiencies( aa5, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( saa5, margEff = TRUE ), 
   efficiencies( saa5i, margEff = TRUE ) )
all.equal( efficiencies( saa5, asInData = TRUE, margEff = TRUE, minusU = FALSE ), 
   efficiencies( saa5i, asInData = TRUE, margEff = TRUE, minusU = FALSE ) )
all.equal( fitted( aa5, asInData = TRUE ) + residuals( aa5, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )
printAll( aa5 )

## cross-section data, efficiency effects frontier, no Z vars
aa9 <- sfa( log( output ) ~ log( capital ) + log( labour ) | - 1,
   data = front41Data )
saa9i <- sfa( logOutput ~ ones + logCapital + logLabour - 1 | - 1, 
   data = front41Data )
all.equal( saa9i[ -c( 3, 7, 20, 42 ) ], aa9[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
round( coef( summary( aa9, extraPar = TRUE ) ), 2 )
print( summary( aa9 ), digits = 1 )
all.equal( summary( aa9 )[ -c( 3, 7, 20, 42 ) ], 
   summary( saa9i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( aa9, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( saa9i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( aa9 )
lrtest( aa9 )
round( efficiencies( aa9 ), 2 )
round( efficiencies( aa9, asInData = TRUE ), 2 )
all.equal( efficiencies( aa9, minusU = FALSE ), 
   efficiencies( saa9i, minusU = FALSE ) )
all.equal( efficiencies( aa9, asInData = TRUE ), 
   efficiencies( saa9i, asInData = TRUE ) )
round( fitted( aa9 ), 2 )
round( fitted( aa9, asInData = TRUE ), 2 )
all.equal( fitted( aa9, asInData = TRUE ) + residuals( aa9, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )


## cross-section data with NAs and infinit values
naData <- front41Data
naData$output[3] <- NA
naData$capital[5] <- 0
naData$labour[9] <- 0
naData$firmNo[14] <- NA

## cross-section data with NAs, error components frontier
San1 <- sfa( log( output ) ~ log( capital ) + log( labour ), data = naData )
San1i <- sfa( log( output ) ~ ones + log( capital ) + log( labour ) - 1, 
   data = naData )
all.equal( San1i[ -c( 3, 7, 20, 42 ) ], San1[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( summary( San1 ), digits = 1 )
all.equal( summary( San1 )[ -c( 3, 7, 20, 42 ) ], 
   summary( San1i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( San1, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( San1i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( San1 )
round( efficiencies( San1 ), 2 )
round( efficiencies( San1, asInData = TRUE ), 2 )
all.equal( efficiencies( San1, minusU = FALSE ), 
   efficiencies( San1i, minusU = FALSE ) )
all.equal( efficiencies( San1, asInData = TRUE ), 
   efficiencies( San1i, asInData = TRUE ) )
round( fitted( San1 ), 2 )
round( fitted( San1, asInData = TRUE ), 2 )
all.equal( fitted( San1, asInData = TRUE ) + residuals( San1, asInData = TRUE ),
   log( naData$output ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( naData$output ) - fitted( San1, asInData = TRUE ),
   residuals( San1, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )

## cross-section data with NAs, efficiency effects frontier
Saan1 <- sfa( log( output ) ~ log( capital ) + log( labour ) | firmNo - 1,
   data = naData )
Saan1i <- sfa( log( output ) ~ ones + log( capital ) + log( labour ) - 1 | 
      firmNo - 1, data = naData )
all.equal( Saan1i[ -c( 3, 7, 20, 42 ) ], Saan1[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( summary( Saan1 ), digits = 1 )
all.equal( summary( Saan1 )[ -c( 3, 7, 20, 42 ) ], 
   summary( Saan1i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( Saan1, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( Saan1i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( Saan1 )
printME( efficiencies( Saan1, margEff = TRUE ) )
printME( efficiencies( Saan1, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( Saan1, margEff = TRUE ), 
   efficiencies( Saan1i, margEff = TRUE ) )
all.equal( efficiencies( Saan1, asInData = TRUE, margEff = TRUE, minusU = FALSE ), 
   efficiencies( Saan1i, asInData = TRUE, margEff = TRUE, minusU = FALSE ) )
round( fitted( Saan1 ), 2 )
round( fitted( Saan1, asInData = TRUE ), 2 )
all.equal( fitted( Saan1, asInData = TRUE ) + residuals( Saan1, asInData = TRUE ),
   log( naData$output ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( naData$output )- fitted( Saan1, asInData = TRUE ),
   residuals( Saan1, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )


## data set of rice producers in the Philippines
data( riceProdPhil )
riceProdPhil$lPROD  <- log( riceProdPhil$PROD )
riceProdPhil$lAREA  <- log( riceProdPhil$AREA )
riceProdPhil$lLABOR <- log( riceProdPhil$LABOR )
riceProdPhil$lNPK   <- log( riceProdPhil$NPK )
riceProdPhil$ones   <- 1

## cross-section rice data, error components frontier
sbb1 <- sfa( lPROD ~ lAREA + lLABOR + lNPK, data = riceProdPhil )
Sbb1 <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ),
   data = riceProdPhil )
all.equal( Sbb1[-42], sbb1[-42], check.attributes = FALSE, tol = 1e-4 )
bb1 <- frontier( data = riceProdPhil,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ) )
all.equal( sbb1[-42], bb1[-42], tol = 1e-4 )
sbb1i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1, data = riceProdPhil )
all.equal( sbb1i[ -c( 3, 7, 20, 42 ) ], sbb1[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( bb1, digits = 1 )
coef( bb1, which = "start" )
round( coef( bb1, which = "ols" ), 2 )
round( coef( bb1, which = "grid" ), 2 )
round( coef( bb1, which = "grid", extraPar = TRUE ), 2 )
round( coef( bb1 ), 2 )
round( coef( bb1, extraPar = TRUE ), 2 )
round( coef( summary( bb1 ), which = "ols" ), 2 )
round( coef( summary( bb1 ) ), 2 )
round( coef( summary( bb1, extraPar = TRUE ) ), 2 )
round( vcov( bb1 ), 2 )
round( vcov( bb1, extraPar = TRUE ), 2 )
nobs( bb1 )
print( summary( bb1 ), digits = 1 )
print( summary( bb1, extraPar = TRUE ), digits = 1 )
all.equal( summary( sbb1, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sbb1i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( sbb1, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sbb1i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( bb1 )
round( efficiencies( bb1 ), 2 )
round( efficiencies( bb1, asInData = TRUE ), 2 )
all.equal( efficiencies( bb1 ), efficiencies( sbb1i ) )
all.equal( efficiencies( bb1, asInData = TRUE, minusU = FALSE ), 
   efficiencies( sbb1i, asInData = TRUE, minusU = FALSE ) )
round( residuals( bb1 ), 2 )
round( residuals( bb1, asInData = TRUE ), 2 )
all.equal( fitted( bb1, asInData = TRUE ) + residuals( bb1, asInData = TRUE ),
   log( riceProdPhil$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( bb1 )

## cross-section rice data, error components frontier, truncNorm
sbb2 <- sfa( lPROD ~ lAREA + lLABOR + lNPK, data = riceProdPhil,
   truncNorm = TRUE )
bb2 <- frontier( data = riceProdPhil,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   truncNorm = TRUE )
all.equal( sbb2[-42], bb2[-42], tol = 1e-4 )
sbb2i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1, data = riceProdPhil,
   truncNorm = TRUE )
all.equal( sbb2i[ -c( 3, 7, 20, 42 ) ], sbb2[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( bb2, digits = 1 )
coef( bb2, which = "start" )
round( coef( bb2, which = "ols" ), 2 )
round( coef( bb2, which = "grid" ), 2 )
round( coef( bb2, which = "grid", extraPar = TRUE ), 2 )
round( coef( bb2 ), 2 )
round( coef( bb2, extraPar = TRUE ), 2 )
round( coef( summary( bb2 ), which = "ols" ), 2 )
round( coef( summary( bb2 ) ), 2 )
round( coef( summary( bb2, extraPar = TRUE ) ), 2 )
round( vcov( bb2 ), 2 )
round( vcov( bb2, extraPar = TRUE ), 2 )
nobs( bb2 )
print( summary( bb2 ), digits = 1 )
print( summary( bb2, effMinusU = FALSE ), digits = 1 )
print( summary( bb2, extraPar = TRUE ), digits = 1 )
all.equal( summary( sbb2, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sbb2i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( sbb2, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sbb2i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( bb2 )
round( efficiencies( bb2 ), 2 )
round( efficiencies( bb2, asInData = TRUE ), 2 )
round( efficiencies( bb2, minusU = FALSE ), 2 )
round( efficiencies( bb2, asInData = TRUE, minusU = FALSE ), 2 )
all.equal( efficiencies( bb2, minusU = FALSE ), 
   efficiencies( sbb2i, minusU = FALSE ) )
all.equal( efficiencies( bb2, asInData = TRUE ), 
   efficiencies( sbb2i, asInData = TRUE ) )
round( residuals( bb2 ), 2 )
round( residuals( bb2, asInData = TRUE ), 2 )
all.equal( fitted( bb2, asInData = TRUE ) + residuals( bb2, asInData = TRUE ),
   log( riceProdPhil$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( bb2 )

## cross-section rice data, efficiency effects frontier
sbb5 <- sfa( lPROD ~ lAREA + lLABOR + lNPK | EDYRS + BANRAT - 1,
   data = riceProdPhil )
Sbb5 <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ) |
   EDYRS + BANRAT - 1, data = riceProdPhil )
all.equal( Sbb5[-42], sbb5[-42], check.attributes = FALSE, tol = 1e-4 )
bb5 <- frontier( data = riceProdPhil,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   zNames = c( "EDYRS", "BANRAT" ) )
all.equal( sbb5[-42], bb5[-42], tol = 1e-4 )
sbb5i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1 | EDYRS + BANRAT - 1, 
   data = riceProdPhil )
all.equal( sbb5i[ -c( 3, 7, 20, 42 ) ], sbb5[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( bb5, digits = 1 )
coef( bb5, which = "start" )
round( coef( bb5, which = "ols" ), 2 )
round( coef( bb5, which = "grid" ), 2 )
round( coef( bb5, which = "grid", extraPar = TRUE ), 2 )
round( coef( bb5 ), 2 )
round( coef( bb5, extraPar = TRUE ), 2 )
round( coef( summary( bb5 ), which = "ols" ), 2 )
round( coef( summary( bb5 ) ), 2 )
round( coef( summary( bb5, extraPar = TRUE ) ), 2 )
round( vcov( bb5 ), 2 )
round( vcov( bb5, extraPar = TRUE ), 2 )
nobs( bb5 )
print( summary( bb5 ), digits = 1 )
print( summary( bb5, extraPar = TRUE ), digits = 1 )
all.equal( summary( sbb5, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sbb5i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( sbb5, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sbb5i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( bb5 )
printME( efficiencies( bb5, margEff = TRUE ) )
printME( efficiencies( bb5, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( bb5 ), efficiencies( sbb5i ) )
all.equal( efficiencies( bb5, asInData = TRUE, minusU = FALSE ), 
   efficiencies( sbb5i, asInData = TRUE, minusU = FALSE ) )
round( residuals( bb5 ), 2 )
round( residuals( bb5, asInData = TRUE ), 2 )
all.equal( fitted( bb5, asInData = TRUE ) + residuals( bb5, asInData = TRUE ),
   log( riceProdPhil$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( bb5 )

## cross-section rice data, efficiency effects frontier, zIntercept
sbb6 <- sfa( lPROD ~ lAREA + lLABOR + lNPK | EDYRS + BANRAT,
   data = riceProdPhil )
bb6 <- frontier( data = riceProdPhil,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   zNames = c( "EDYRS", "BANRAT" ), zIntercept = TRUE )
all.equal( sbb6[-42], bb6[-42], tol = 1e-4 )
sbb6i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1 | EDYRS + BANRAT, 
   data = riceProdPhil )
all.equal( sbb6i[ -c( 3, 7, 20, 42 ) ], sbb6[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( bb6, digits = 1 )
coef( bb6, which = "start" )
round( coef( bb6, which = "ols" ), 2 )
round( coef( bb6, which = "grid" ), 2 )
round( coef( bb6, which = "grid", extraPar = TRUE ), 2 )
round( coef( bb6 ), 2 )
round( coef( bb6, extraPar = TRUE ), 2 )
round( coef( summary( bb6 ), which = "ols" ), 2 )
round( coef( summary( bb6 ) ), 2 )
round( coef( summary( bb6, extraPar = TRUE ) ), 2 )
round( vcov( bb6 ), 2 )
round( vcov( bb6, extraPar = TRUE ), 2 )
nobs( bb6 )
print( summary( bb6 ), digits = 1 )
print( summary( bb6, effMinusU = FALSE ), digits = 1 )
print( summary( bb6, extraPar = TRUE ), digits = 1 )
all.equal( summary( sbb6, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sbb6i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ],
   check.attributes = FALSE )
all.equal( summary( sbb6, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sbb6i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( bb6 )
printME( bb6eff <- efficiencies( bb6, margEff = TRUE ) )
printME( bb6effD <- efficiencies( bb6, asInData = TRUE, margEff = TRUE ) )
printME( bb6effF <- efficiencies( bb6, minusU = FALSE, margEff = TRUE ) )
printME( bb6effDF <- efficiencies( bb6, asInData = TRUE, minusU = FALSE, 
   margEff = TRUE ) )
bb6m1 <- bb6
bb6m1$dataTable[ , "EDYRS" ] <- bb6m1$dataTable[ , "EDYRS" ] + 1e-6
all.equal( attr( bb6eff, "margEff" )[ , 1, 1 ], 
   ( efficiencies( bb6m1 ) - bb6eff )[ , 1 ] / 1e-6, tol = 1e-4 )
all.equal( attr( bb6effD, "margEff" )[ , 1 ], 
   c( efficiencies( bb6m1, asInData = TRUE ) - bb6effD ) / 1e-6, tol = 1e-4 )
all.equal( attr( bb6effF, "margEff" )[ , 1, 1 ], 
   ( efficiencies( bb6m1, minusU = FALSE ) - bb6effF )[ , 1 ] / 1e-6, tol = 1e-4 )
all.equal( attr( bb6effDF, "margEff" )[ , 1 ],
   c( efficiencies( bb6m1, asInData = TRUE, minusU = FALSE ) - bb6effDF ) / 1e-6,
   tol = 1e-4 )
bb6m2 <- bb6
bb6m2$dataTable[ , "BANRAT" ] <- bb6m2$dataTable[ , "BANRAT" ] + 1e-6
all.equal( attr( bb6eff, "margEff" )[ , 1, 2 ], 
   ( efficiencies( bb6m2 ) - bb6eff )[ , 1 ] / 1e-6, tol = 1e-4 )
all.equal( attr( bb6effD, "margEff" )[ , 2 ], 
   c( efficiencies( bb6m2, asInData = TRUE ) - bb6effD ) / 1e-6, tol = 1e-4 )
all.equal( attr( bb6effF, "margEff" )[ , 1, 2 ], 
   ( efficiencies( bb6m2, minusU = FALSE ) - bb6effF )[ , 1 ] / 1e-6, tol = 1e-4 )
all.equal( attr( bb6effDF, "margEff" )[ , 2 ],
   c( efficiencies( bb6m2, asInData = TRUE, minusU = FALSE ) - bb6effDF ) / 1e-6,
   tol = 1e-4 )
all.equal( efficiencies( bb6, margEff = TRUE, minusU = FALSE ), 
   efficiencies( sbb6i, margEff = TRUE, minusU = FALSE ) )
all.equal( efficiencies( bb6, asInData = TRUE, margEff = TRUE ), 
   efficiencies( sbb6i, asInData = TRUE, margEff = TRUE ) )
round( residuals( bb6 ), 2 )
round( residuals( bb6, asInData = TRUE ), 2 )
all.equal( fitted( bb6, asInData = TRUE ) + residuals( bb6, asInData = TRUE ),
   log( riceProdPhil$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( bb6 )

## cross-section rice data, error components frontier, truncNorm, starting values
sbb7 <- sfa( lPROD ~ lAREA + lLABOR + lNPK, data = riceProdPhil,
   truncNorm = TRUE, startVal = c( -1, 0.3, 0.3, 0.3, 0.2, 0.9, -0.01 ) )
bb7 <- frontier( data = riceProdPhil,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   truncNorm = TRUE,
   startVal = c( -1, 0.3, 0.3, 0.3, 0.2, 0.9, -0.01 ) )
all.equal( sbb7[-42], bb7[-42], tol = 1e-4 )
sbb7i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1, data = riceProdPhil, 
   truncNorm = TRUE, startVal = c( -1, 0.3, 0.3, 0.3, 0.2, 0.9, -0.01 ) )
all.equal( sbb7i[ -c( 3, 7, 21, 42 ) ], sbb7[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( bb7, digits = 1 )
coef( bb7, which = "start" )
round( coef( bb7, which = "ols" ), 2 )
coef( bb7, which = "grid" )
coef( bb7, which = "grid", extraPar = TRUE )
round( coef( bb7 ), 2 )
round( coef( bb7, extraPar = TRUE ), 2 )
round( coef( summary( bb7 ), which = "ols" ), 2 )
round( coef( summary( bb7 ) ), 2 )
round( coef( summary( bb7, extraPar = TRUE ) ), 2 )
round( vcov( bb7 ), 2 )
round( vcov( bb7, extraPar = TRUE ), 2 )
nobs( bb7 )
print( summary( bb7 ), digits = 1 )
print( summary( bb7, extraPar = TRUE ), digits = 1 )
all.equal( summary( sbb7, extraPar = TRUE )[ -c( 3, 7, 21, 42 ) ], 
   summary( sbb7i, extraPar = TRUE )[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( sbb7, effMinusU = FALSE )[ -c( 3, 7, 21, 42 ) ], 
   summary( sbb7i, effMinusU = FALSE )[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE )
lrtest( bb7 )
round( efficiencies( bb7 ), 2 )
round( efficiencies( bb7, asInData = TRUE ), 2 )
all.equal( efficiencies( bb7 ), efficiencies( sbb7i ) )
all.equal( efficiencies( bb7, asInData = TRUE, minusU = FALSE ), 
   efficiencies( sbb7i, asInData = TRUE, minusU = FALSE ) )
all.equal( fitted( bb7, asInData = TRUE ) + residuals( bb7, asInData = TRUE ),
   log( riceProdPhil$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( bb7 )

## cross-section rice data, efficiency effects frontier, zIntercept, starting values
sbb8 <- sfa( lPROD ~ lAREA + lLABOR + lNPK | EDYRS + BANRAT,
   data = riceProdPhil,
   startVal = c( -1, 0.3, 0.3, 0.3, -0.2, -0.01, -0.3, 0.3, 0.8 ) )
bb8 <- frontier( data = riceProdPhil,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   zNames = c( "EDYRS", "BANRAT" ), zIntercept = TRUE,
   startVal = c( -1, 0.3, 0.3, 0.3, -0.2, -0.01, -0.3, 0.3, 0.8 ) )
all.equal( sbb8[-42], bb8[-42], tol = 1e-4 )
sbb8i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1 | EDYRS + BANRAT, 
   data = riceProdPhil, 
   startVal = c( -1, 0.3, 0.3, 0.3, -0.2, -0.01, -0.3, 0.3, 0.8 ) )
all.equal( sbb8i[ -c( 3, 7, 21, 42 ) ], sbb8[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( bb8, digits = 1 )
coef( bb8, which = "start" )
round( coef( bb8, which = "ols" ), 2 )
coef( bb8, which = "grid" )
coef( bb8, which = "grid", extraPar = TRUE )
round( coef( bb8 ), 2 )
round( coef( bb8, extraPar = TRUE ), 2 )
round( coef( summary( bb8 ), which = "ols" ), 2 )
round( coef( summary( bb8 ) ), 2 )
round( coef( summary( bb8, extraPar = TRUE ) ), 2 )
round( vcov( bb8 ), 2 )
round( vcov( bb8, extraPar = TRUE ), 2 )
nobs( bb8 )
print( summary( bb8 ), digits = 1 )
print( summary( bb8, extraPar = TRUE ), digits = 1 )
all.equal( summary( sbb8, extraPar = TRUE )[ -c( 3, 7, 21, 42 ) ], 
   summary( sbb8i, extraPar = TRUE )[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( sbb8, effMinusU = FALSE )[ -c( 3, 7, 21, 42 ) ], 
   summary( sbb8i, effMinusU = FALSE )[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE )
lrtest( bb8 )
printME( efficiencies( bb8, margEff = TRUE ) )
printME( efficiencies( bb8, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( bb8, margEff = TRUE, minusU = FALSE ), 
   efficiencies( sbb8i, margEff = TRUE, minusU = FALSE ) )
all.equal( efficiencies( bb8, asInData = TRUE, margEff = TRUE ), 
   efficiencies( sbb8i, asInData = TRUE, margEff = TRUE ) )
all.equal( fitted( bb8, asInData = TRUE ) + residuals( bb8, asInData = TRUE ),
   log( riceProdPhil$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( bb8 )

## cross-section rice data, efficiency effects frontier, no Z vars
bb9 <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ) | - 1,
   data = riceProdPhil )
sbb9i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1 | - 1, 
   data = riceProdPhil )
all.equal( sbb9i[ -c( 3, 7, 20, 42 ) ], bb9[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
round( coef( summary( bb9, extraPar = TRUE ) ), 2 )
print( summary( bb9 ), digits = 1 )
all.equal( summary( bb9 )[ -c( 3, 7, 20, 42 ) ], 
   summary( sbb9i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( bb9, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sbb9i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( bb9 )
lrtest( bb9 )
all.equal( efficiencies( bb9, minusU = FALSE ), 
   efficiencies( sbb9i, minusU = FALSE ) )
all.equal( efficiencies( bb9, asInData = TRUE ), 
   efficiencies( sbb9i, asInData = TRUE ) )
round( residuals( bb9 ), 2 )
round( residuals( bb9, asInData = TRUE ), 2 )
all.equal( fitted( bb9, asInData = TRUE ) + residuals( bb9, asInData = TRUE ),
   log( riceProdPhil$PROD ), check.attributes = FALSE, tol = 1e-4 )


## Cost Frontier (with land as quasi-fixed input)
riceProdPhil$cost <- riceProdPhil$LABOR * riceProdPhil$LABORP +
   riceProdPhil$NPK * riceProdPhil$NPKP
riceProdPhil$lCost   <- log( riceProdPhil$cost )
riceProdPhil$lLABORP <- log( riceProdPhil$LABORP )
riceProdPhil$lNPKP   <- log( riceProdPhil$NPKP )

## cross-section rice data, error components cost frontier
sdd1 <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP, data = riceProdPhil,
   ineffDecrease = FALSE )
Sdd1 <- sfa( log( cost ) ~ log( PROD ) + log( AREA ) + log( LABORP ) +
   log( NPKP ), data = riceProdPhil, ineffDecrease = FALSE )
all.equal( Sdd1[-42], sdd1[-42], check.attributes = FALSE, tol = 1e-4 )
dd1 <- frontier( "lCost", xNames = c( "lPROD", "lAREA", "lLABORP", "lNPKP" ),
   data = riceProdPhil, ineffDecrease = FALSE )
all.equal( sdd1[-42], dd1[-42], tol = 1e-4 )
sdd1i <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1, 
   data = riceProdPhil, ineffDecrease = FALSE )
all.equal( sdd1i[ -c( 3, 7, 20, 42 ) ], sdd1[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( dd1, digits = 1 )
coef( dd1, which = "start" )
round( coef( dd1, which = "ols" ), 2 )
round( coef( dd1, which = "grid" ), 2 )
round( coef( dd1, which = "grid", extraPar = TRUE ), 2 )
round( coef( dd1 ), 2 )
round( coef( dd1, extraPar = TRUE ), 2 )
round( coef( summary( dd1 ), which = "ols" ), 2 )
round( coef( summary( dd1 ) ), 2 )
round( coef( summary( dd1, extraPar = TRUE ) ), 2 )
round( vcov( dd1 ), 2 )
round( vcov( dd1, extraPar = TRUE ), 2 )
nobs( dd1 )
print( summary( dd1 ), digits = 1 )
print( summary( dd1, effMinusU = FALSE ), digits = 1 )
print( summary( dd1, extraPar = TRUE ), digits = 1 )
all.equal( summary( sdd1, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sdd1i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( sdd1, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sdd1i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( dd1 )
round( efficiencies( dd1 ), 2 )
round( efficiencies( dd1, asInData = TRUE ), 2 )
round( efficiencies( dd1, minusU = FALSE ), 2 )
round( efficiencies( dd1, asInData = TRUE, minusU = FALSE ), 2 )
all.equal( efficiencies( dd1 ), efficiencies( sdd1i ) )
all.equal( efficiencies( dd1, asInData = TRUE, minusU = FALSE ), 
   efficiencies( sdd1i, asInData = TRUE, minusU = FALSE ) )
round( residuals( dd1 ), 2 )
round( residuals( dd1, asInData = TRUE ), 2 )
all.equal( fitted( dd1, asInData = TRUE ) + residuals( dd1, asInData = TRUE ),
   log( riceProdPhil$cost ), check.attributes = FALSE, tol = 1e-4 )
printAll( dd1 )

## cross-section rice data, error components cost frontier, truncNorm
sdd2 <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP, data = riceProdPhil,
   ineffDecrease = FALSE, truncNorm = TRUE )
dd2 <- frontier( "lCost", xNames = c( "lPROD", "lAREA", "lLABORP", "lNPKP" ),
   data = riceProdPhil, ineffDecrease = FALSE, truncNorm = TRUE )
all.equal( sdd2[-42], dd2[-42], tol = 1e-4 )
sdd2i <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1, 
   data = riceProdPhil, ineffDecrease = FALSE, truncNorm = TRUE )
all.equal( sdd2i[ -c( 3, 7, 20, 42 ) ], sdd2[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( dd2, digits = 1 )
coef( dd2, which = "start" )
round( coef( dd2, which = "ols" ), 2 )
round( coef( dd2, which = "grid" ), 2 )
round( coef( dd2, which = "grid", extraPar = TRUE ), 2 )
round( coef( dd2 ), 2 )
round( coef( dd2, extraPar = TRUE ), 2 )
round( coef( summary( dd2 ), which = "ols" ), 2 )
round( coef( summary( dd2 ) ), 2 )
round( coef( summary( dd2, extraPar = TRUE ) ), 2 )
round( vcov( dd2 ), 2 )
round( vcov( dd2, extraPar = TRUE ), 2 )
nobs( dd2 )
print( summary( dd2, effMinusU = FALSE ), digits = 1 )
print( summary( dd2, extraPar = TRUE ), digits = 1 )
all.equal( summary( sdd2, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sdd2i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( sdd2, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sdd2i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( dd2 )
round( efficiencies( dd2, minusU = FALSE ), 2 )
round( efficiencies( dd2, asInData = TRUE, minusU = FALSE ), 2 )
all.equal( efficiencies( dd2 ), efficiencies( sdd2i ) )
all.equal( efficiencies( dd2, asInData = TRUE, minusU = FALSE ), 
   efficiencies( sdd2i, asInData = TRUE, minusU = FALSE ) )
round( residuals( dd2 ), 2 )
round( residuals( dd2, asInData = TRUE ), 2 )
all.equal( fitted( dd2, asInData = TRUE ) + residuals( dd2, asInData = TRUE ),
   log( riceProdPhil$cost ), check.attributes = FALSE, tol = 1e-4 )
printAll( dd2 )

## cross-section rice data, efficiency effects cost frontier
sdd5 <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP | EDYRS + BANRAT - 1,
   data = riceProdPhil, ineffDecrease = FALSE )
Sdd5 <- sfa( log( cost ) ~ log( PROD ) + log( AREA ) + log( LABORP ) +
   log( NPKP ) | EDYRS + BANRAT - 1, data = riceProdPhil, ineffDecrease = FALSE )
all.equal( Sdd5[-42], sdd5[-42], check.attributes = FALSE, tol = 1e-4 )
dd5 <- frontier( "lCost", xNames = c( "lPROD", "lAREA", "lLABORP", "lNPKP" ),
   zNames = c( "EDYRS", "BANRAT" ), data = riceProdPhil,
   ineffDecrease = FALSE )
all.equal( sdd5[-42], dd5[-42], tol = 1e-4 )
sdd5i <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1 |
      EDYRS + BANRAT - 1, data = riceProdPhil, ineffDecrease = FALSE )
all.equal( sdd5i[ -c( 3, 7, 20, 42 ) ], sdd5[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( dd5, digits = 1 )
coef( dd5, which = "start" )
round( coef( dd5, which = "ols" ), 2 )
round( coef( dd5, which = "grid" ), 2 )
round( coef( dd5, which = "grid", extraPar = TRUE ), 2 )
round( coef( dd5 ), 2 )
round( coef( dd5, extraPar = TRUE ), 2 )
round( coef( summary( dd5 ), which = "ols" ), 2 )
round( coef( summary( dd5 ) ), 2 )
round( coef( summary( dd5, extraPar = TRUE ) ), 2 )
round( vcov( dd5 ), 2 )
round( vcov( dd5, extraPar = TRUE ), 2 )
nobs( dd5 )
print( summary( dd5 ), digits = 1 )
print( summary( dd5, effMinusU = FALSE ), digits = 1 )
print( summary( dd5, extraPar = TRUE ), digits = 1 )
all.equal( summary( sdd5, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sdd5i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( sdd5, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sdd5i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( dd5 )
printME( dd5eff <- efficiencies( dd5, margEff = TRUE ) )
printME( dd5effD <- efficiencies( dd5, asInData = TRUE, margEff = TRUE ) )
printME( dd5effF <- efficiencies( dd5, minusU = FALSE, margEff = TRUE ) )
printME( dd5effDF <- efficiencies( dd5, asInData = TRUE, minusU = FALSE, 
   margEff = TRUE ) )
dd5m1 <- dd5
dd5m1$dataTable[ , "EDYRS" ] <- dd5m1$dataTable[ , "EDYRS" ] + 1e-6
all.equal( attr( dd5eff, "margEff" )[ , 1, 1 ], 
   ( efficiencies( dd5m1 ) - dd5eff )[ , 1 ] / 1e-6, tol = 1e-4 )
all.equal( attr( dd5effD, "margEff" )[ , 1 ], 
   c( efficiencies( dd5m1, asInData = TRUE ) - dd5effD ) / 1e-6, tol = 1e-4 )
all.equal( attr( dd5effF, "margEff" )[ , 1, 1 ], 
   ( efficiencies( dd5m1, minusU = FALSE ) - dd5effF )[ , 1 ] / 1e-6, tol = 1e-4 )
all.equal( attr( dd5effDF, "margEff" )[ , 1 ],
   c( efficiencies( dd5m1, asInData = TRUE, minusU = FALSE ) - dd5effDF ) / 1e-6,
   tol = 1e-4 )
dd5m2 <- dd5
dd5m2$dataTable[ , "BANRAT" ] <- dd5m2$dataTable[ , "BANRAT" ] + 1e-6
all.equal( attr( dd5eff, "margEff" )[ , 1, 2 ], 
   ( efficiencies( dd5m2 ) - dd5eff )[ , 1 ] / 1e-6, tol = 1e-4 )
all.equal( attr( dd5effD, "margEff" )[ , 2 ], 
   c( efficiencies( dd5m2, asInData = TRUE ) - dd5effD ) / 1e-6, tol = 1e-4 )
all.equal( attr( dd5effF, "margEff" )[ , 1, 2 ], 
   ( efficiencies( dd5m2, minusU = FALSE ) - dd5effF )[ , 1 ] / 1e-6, tol = 1e-4 )
all.equal( attr( dd5effDF, "margEff" )[ , 2 ],
   c( efficiencies( dd5m2, asInData = TRUE, minusU = FALSE ) - dd5effDF ) / 1e-6,
   tol = 1e-4 )
all.equal( efficiencies( dd5, margEff = TRUE, minusU = FALSE ), 
   efficiencies( sdd5i, margEff = TRUE, minusU = FALSE ) )
all.equal( efficiencies( dd5, asInData = TRUE, margEff = TRUE ), 
   efficiencies( sdd5i, asInData = TRUE, margEff = TRUE ) )
round( residuals( dd5 ), 2 )
round( residuals( dd5, asInData = TRUE ), 2 )
all.equal( fitted( dd5, asInData = TRUE ) + residuals( dd5, asInData = TRUE ),
   log( riceProdPhil$cost ), check.attributes = FALSE, tol = 1e-4 )
printAll( dd5 )

## cross-section rice data, efficiency effects cost frontier, zIntercept
sdd6 <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP | EDYRS + BANRAT,
   data = riceProdPhil, ineffDecrease = FALSE )
dd6 <- frontier( "lCost", xNames = c( "lPROD", "lAREA", "lLABORP", "lNPKP" ),
   zNames = c( "EDYRS", "BANRAT" ), data = riceProdPhil,
   ineffDecrease = FALSE, zIntercept = TRUE )
all.equal( sdd6[-42], dd6[-42], tol = 1e-4 )
sdd6i <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1 |
      EDYRS + BANRAT, data = riceProdPhil, ineffDecrease = FALSE )
all.equal( sdd6i[ -c( 3, 7, 20, 42 ) ], sdd6[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( dd6, digits = 1 )
coef( dd6, which = "start" )
round( coef( dd6, which = "ols" ), 2 )
round( coef( dd6, which = "grid" ), 2 )
round( coef( dd6, which = "grid", extraPar = TRUE ), 2 )
round( coef( dd6 ), 2 )
round( coef( dd6, extraPar = TRUE ), 2 )
round( coef( summary( dd6 ), which = "ols" ), 2 )
round( coef( summary( dd6 ) ), 2 )
round( coef( summary( dd6, extraPar = TRUE ) ), 2 )
round( vcov( dd6 ), 2 )
round( vcov( dd6, extraPar = TRUE ), 2 )
nobs( dd6 )
print( summary( dd6, effMinusU = FALSE ), digits = 1 )
print( summary( dd6, extraPar = TRUE ), digits = 1 )
all.equal( summary( sdd6, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sdd6i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( sdd6, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sdd6i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( dd6 )
printME( efficiencies( dd6, minusU = FALSE, margEff = TRUE ) )
printME( efficiencies( dd6, asInData = TRUE, minusU = FALSE, margEff = TRUE ) )
all.equal( efficiencies( dd6, margEff = TRUE ), 
   efficiencies( sdd6i, margEff = TRUE ) )
all.equal( efficiencies( dd6, asInData = TRUE, margEff = TRUE, minusU = FALSE ), 
   efficiencies( sdd6i, asInData = TRUE, margEff = TRUE, minusU = FALSE ) )
round( residuals( dd6 ), 2 )
round( residuals( dd6, asInData = TRUE ), 2 )
all.equal( fitted( dd6, asInData = TRUE ) + residuals( dd6, asInData = TRUE ),
   log( riceProdPhil$cost ), check.attributes = FALSE, tol = 1e-4 )
printAll( dd6 )

## cross-section rice data, efficiency effects cost frontier: no Z vars
dd9 <- sfa( log( cost ) ~ log( PROD ) + log( AREA ) + log( LABORP ) +
   log( NPKP ) | - 1, data = riceProdPhil, ineffDecrease = FALSE )
sdd9i <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1 | - 1, 
   data = riceProdPhil, ineffDecrease = FALSE )
all.equal( sdd9i[ -c( 3, 7, 20, 42 ) ], dd9[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
round( coef( summary( dd9, extraPar = TRUE ) ), 2 )
print( summary( dd9, effMinusU = FALSE ), digits = 1 )
all.equal( summary( dd9 )[ -c( 3, 7, 20, 42 ) ], 
   summary( sdd9i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( dd9, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sdd9i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( dd9 )
lrtest( dd9 )
round( efficiencies( dd9 ), 2 )
round( efficiencies( dd9, asInData = TRUE ), 2 )
all.equal( efficiencies( dd9, minusU = FALSE ), 
   efficiencies( sdd9i, minusU = FALSE ) )
all.equal( efficiencies( dd9, asInData = TRUE ), 
   efficiencies( sdd9i, asInData = TRUE ) )
round( fitted( dd9 ), 2 )
round( fitted( dd9, asInData = TRUE ), 2 )
all.equal( fitted( dd9, asInData = TRUE ) + residuals( dd9, asInData = TRUE ),
   log( riceProdPhil$cost ), check.attributes = FALSE, tol = 1e-4 )


## error components frontier
## with "true" fixed individual effects and observation-specific efficiencies
riceTrue <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ) + 
   factor( FMERCODE ),  data = riceProdPhil )
summary( riceTrue )
lrtest( riceTrue )
lrtest( Sbb1, riceTrue )
round( fitted( riceTrue ), 2 )
round( fitted( riceTrue, asInData = TRUE ), 2 )
all.equal( fitted( riceTrue, asInData = TRUE ) + residuals( riceTrue, asInData = TRUE ),
   log( riceProdPhil$PROD ), check.attributes = FALSE, tol = 1e-4 )


## panel data
riceProdPhil$farm <- paste( "F_", ifelse( riceProdPhil$FMERCODE > 9, "", "0" ),
   riceProdPhil$FMERCODE, sep = "" )
riceProdPhil$year <- riceProdPhil$YEARDUM + 1998
riceProdPhilPanel <- plm.data( riceProdPhil, c( "farm", "year" ) )

## panel data, error components frontier
sb1 <- sfa( lPROD ~ lAREA + lLABOR + lNPK, data = riceProdPhilPanel,
   printIter = 2 )
Sb1 <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ),
   data = riceProdPhilPanel )
all.equal( Sb1[-c(8,42)], sb1[-c(8,42)], check.attributes = FALSE, tol = 1e-4 )
b1 <- frontier( data = riceProdPhilPanel,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ) )
all.equal( sb1[-c(8,42)], b1[-c(8,42)], tol = 1e-4 )
sb1i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1, 
   data = riceProdPhilPanel )
all.equal( sb1i[ -c( 3, 7, 8, 20, 42 ) ], sb1[ -c( 3, 7, 8, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b1, digits = 1 )
coef( b1, which = "start" )
round( coef( b1, which = "ols" ), 2 )
round( coef( b1, which = "grid" ), 2 )
round( coef( b1, which = "grid", extraPar = TRUE ), 2 )
round( coef( b1 ), 2 )
round( coef( b1, extraPar = TRUE ), 2 )
round( coef( summary( b1 ), which = "ols" ), 2 )
round( coef( summary( b1 ) ), 2 )
round( coef( summary( b1, extraPar = TRUE ) ), 2 )
round( vcov( b1 ), 2 )
round( vcov( b1, extraPar = TRUE ), 2 )
print( logLik( b1, which = "ols" ), digits = 4 )
print( logLik( b1 ), digits = 4 )
nobs( b1 )
print( summary( b1 ), digits = 1 )
print( summary( b1, effMinusU = FALSE ), digits = 1 )
print( summary( b1, extraPar = TRUE ), digits = 1 )
all.equal( summary( b1, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb1i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ],
   check.attributes = FALSE )
all.equal( summary( b1, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb1i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( b1 )
round( efficiencies( b1 ), 2 )
round( efficiencies( b1, asInData = TRUE ), 2 )
round( efficiencies( b1, minusU = FALSE ), 2 )
round( efficiencies( b1, asInData = TRUE, minusU = FALSE ), 2 )
all.equal( efficiencies( b1 ), efficiencies( sb1i ) )
all.equal( efficiencies( b1, asInData = TRUE, minusU = FALSE ), 
   efficiencies( sb1i, asInData = TRUE, minusU = FALSE ) )
round( residuals( b1 ), 2 )
round( residuals( b1, asInData = TRUE ), 2 )
all.equal( fitted( b1, asInData = TRUE ) + residuals( b1, asInData = TRUE ),
   log( riceProdPhilPanel$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( b1 )

## panel data, error components frontier, truncNorm
sb2 <- sfa( lPROD ~ lAREA + lLABOR + lNPK, data = riceProdPhilPanel,
   truncNorm = TRUE )
b2 <- frontier( data = riceProdPhilPanel,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   truncNorm = TRUE )
all.equal( sb2[-42], b2[-42], tol = 1e-4 )
sb2i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1, 
   data = riceProdPhilPanel, truncNorm = TRUE )
all.equal( sb2i[ -c( 3, 7, 20, 42 ) ], sb2[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b2, digits = 1 )
coef( b2, which = "start" )
round( coef( b2, which = "ols" ), 2 )
round( coef( b2, which = "grid" ), 2 )
round( coef( b2, which = "grid", extraPar = TRUE ), 2 )
round( coef( b2 ), 2 )
round( coef( b2, extraPar = TRUE ), 2 )
round( coef( summary( b2 ), which = "ols" ), 2 )
round( coef( summary( b2 ) ), 2 )
round( coef( summary( b2, extraPar = TRUE ) ), 2 )
round( vcov( b2 ), 2 )
round( vcov( b2, extraPar = TRUE ), 2 )
print( logLik( b2, which = "ols" ), digits = 4 )
print( logLik( b2 ), digits = 4 )
nobs( b2 )
print( summary( b2 ), digits = 1 )
print( summary( b2, extraPar = TRUE ), digits = 1 )
all.equal( summary( b2, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb2i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( b2, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb2i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( b2 )
round( efficiencies( b2 ), 2 )
round( efficiencies( b2, asInData = TRUE ), 2 )
all.equal( efficiencies( b2, minusU = FALSE ), 
   efficiencies( sb2i, minusU = FALSE ) )
all.equal( efficiencies( b2, asInData = TRUE ), 
   efficiencies( sb2i, asInData = TRUE ) )
round( fitted( b2 ), 2 )
round( fitted( b2, asInData = TRUE ), 2 )
round( residuals( b2 ), 2 )
round( residuals( b2, asInData = TRUE ), 2 )
all.equal( fitted( b2, asInData = TRUE ) + residuals( b2, asInData = TRUE ),
   log( riceProdPhilPanel$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( b2 )

## panel data, error components frontier, timeEffect
sb3 <- sfa( lPROD ~ lAREA + lLABOR + lNPK, data = riceProdPhilPanel,
   timeEffect = TRUE )
b3 <- frontier( data = riceProdPhilPanel,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   timeEffect = TRUE )
all.equal( sb3[-42], b3[-42], tol = 1e-4 )
sb3i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1, 
   data = riceProdPhilPanel, timeEffect = TRUE )
all.equal( sb3i[ -c( 3, 7, 20, 42 ) ], sb3[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b3, digits = 1 )
coef( b3, which = "start" )
round( coef( b3, which = "ols" ), 2 )
round( coef( b3, which = "grid" ), 2 )
round( coef( b3, which = "grid", extraPar = TRUE ), 2 )
round( coef( b3 ), 2 )
round( coef( b3, extraPar = TRUE ), 2 )
round( coef( summary( b3 ), which = "ols" ), 2 )
round( coef( summary( b3 ) ), 2 )
round( coef( summary( b3, extraPar = TRUE ) ), 2 )
round( vcov( b3 ), 2 )
round( vcov( b3, extraPar = TRUE ), 2 )
print( logLik( b3, which = "ols" ), digits = 4 )
print( logLik( b3 ), digits = 4 )
nobs( b3 )
print( summary( b3 ), digits = 1 )
print( summary( b3, extraPar = TRUE ), digits = 1 )
all.equal( summary( b3, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb3i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( b3, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb3i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( b3 )
round( efficiencies( b3 ), 2 )
round( efficiencies( b3, asInData = TRUE ), 2 )
all.equal( efficiencies( b3 ), efficiencies( sb3i ) )
all.equal( efficiencies( b3, asInData = TRUE, minusU = FALSE ), 
   efficiencies( sb3i, asInData = TRUE, minusU = FALSE ) )
round( residuals( b3 ), 2 )
round( residuals( b3, asInData = TRUE ), 2 )
all.equal( fitted( b3, asInData = TRUE ) + residuals( b3, asInData = TRUE ),
   log( riceProdPhilPanel$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( b3 )

## panel data, error components frontier, truncNorm, timeEffect
sb4 <- sfa( lPROD ~ lAREA + lLABOR + lNPK, data = riceProdPhilPanel,
   truncNorm = TRUE, timeEffect = TRUE )
b4 <- frontier( data = riceProdPhilPanel,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   truncNorm = TRUE, timeEffect = TRUE )
all.equal( sb4[-42], b4[-42], tol = 1e-4 )
sb4i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1, 
   data = riceProdPhilPanel, truncNorm = TRUE, timeEffect = TRUE )
all.equal( sb4i[ -c( 3, 7, 20, 42 ) ], sb4[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b4, digits = 1 )
coef( b4, which = "start" )
round( coef( b4, which = "ols" ), 2 )
round( coef( b4, which = "grid" ), 2 )
round( coef( b4, which = "grid", extraPar = TRUE ), 2 )
round( coef( b4 ), 2 )
round( coef( b4, extraPar = TRUE ), 2 )
round( coef( summary( b4 ), which = "ols" ), 2 )
round( coef( summary( b4 ) ), 2 )
round( coef( summary( b4, extraPar = TRUE ) ), 2 )
round( vcov( b4 ), 2 )
round( vcov( b4, extraPar = TRUE ), 2 )
print( logLik( b4, which = "ols" ), digits = 4 )
print( logLik( b4 ), digits = 4 )
nobs( b4 )
print( summary( b4 ), digits = 1 )
print( summary( b4, extraPar = TRUE ), digits = 1 )
all.equal( summary( b4, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb4i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( b4, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb4i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( b4 )
round( efficiencies( b4 ), 2 )
round( efficiencies( b4, asInData = TRUE ), 2 )
all.equal( efficiencies( b4, minusU = FALSE ), 
   efficiencies( sb4i, minusU = FALSE ) )
all.equal( efficiencies( b4, asInData = TRUE ), 
   efficiencies( sb4i, asInData = TRUE ) )
round( residuals( b4 ), 2 )
round( residuals( b4, asInData = TRUE ), 2 )
all.equal( fitted( b4, asInData = TRUE ) + residuals( b4, asInData = TRUE ),
   log( riceProdPhilPanel$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( b4 )

## error components frontier
## with "true" fixed individual effects and time-invariant efficiencies
ricePanelTrue <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ) + 
      factor( FMERCODE ),  data = riceProdPhilPanel )
summary( ricePanelTrue )
lrtest( ricePanelTrue )
lrtest( Sb1, ricePanelTrue )
round( fitted( ricePanelTrue ), 2 )
round( fitted( ricePanelTrue, asInData = TRUE ), 2 )
all.equal( fitted( ricePanelTrue, asInData = TRUE ) + 
      residuals( ricePanelTrue, asInData = TRUE ),
   log( riceProdPhilPanel$PROD ), check.attributes = FALSE, tol = 1e-4 )

## error components frontier
## with "true" fixed individual effects and time-variant efficiencies
ricePanelTimeTrue <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ) + 
      factor( FMERCODE ), data = riceProdPhilPanel, timeEffect = TRUE )
summary( ricePanelTimeTrue )
lrtest( ricePanelTimeTrue )
lrtest( sb3, ricePanelTimeTrue )
lrtest( ricePanelTrue, ricePanelTimeTrue )
all.equal( fitted( ricePanelTimeTrue, asInData = TRUE ) + 
      residuals( ricePanelTimeTrue, asInData = TRUE ),
   log( riceProdPhilPanel$PROD ), check.attributes = FALSE, tol = 1e-4 )


## panel data, efficiency effects frontier
sb5 <- sfa( lPROD ~ lAREA + lLABOR + lNPK | EDYRS + BANRAT - 1,
   data = riceProdPhilPanel )
Sb5 <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ) |
   EDYRS + BANRAT - 1, data = riceProdPhilPanel, printIter = 5 )
all.equal( Sb5[-c(8,42)], sb5[-c(8,42)], check.attributes = FALSE, tol = 1e-4 )
b5 <- frontier( data = riceProdPhilPanel,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   zNames = c( "EDYRS", "BANRAT" ) )
all.equal( sb5[-42], b5[-42], tol = 1e-4 )
sb5i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1 | EDYRS + BANRAT - 1, 
   data = riceProdPhilPanel )
all.equal( sb5i[ -c( 3, 7, 20, 42 ) ], sb5[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b5, digits = 1 )
coef( b5, which = "start" )
round( coef( b5, which = "ols" ), 2 )
round( coef( b5, which = "grid" ), 2 )
round( coef( b5, which = "grid", extraPar = TRUE ), 2 )
round( coef( b5 ), 2 )
round( coef( b5, extraPar = TRUE ), 2 )
round( coef( summary( b5 ), which = "ols" ), 2 )
round( coef( summary( b5 ) ), 2 )
round( coef( summary( b5, extraPar = TRUE ) ), 2 )
round( vcov( b5 ), 2 )
round( vcov( b5, extraPar = TRUE ), 2 )
print( logLik( b5, which = "ols" ), digits = 4 )
print( logLik( b5 ), digits = 4 )
nobs( b5 )
print( summary( b5 ), digits = 1 )
print( summary( b5, effMinusU = FALSE ), digits = 1 )
print( summary( b5, extraPar = TRUE ), digits = 1 )
all.equal( summary( b5, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb5i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( b5, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb5i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( b5 )
printME( b5eff <- efficiencies( b5, margEff = TRUE ) )
printME( b5effD <- efficiencies( b5, asInData = TRUE, margEff = TRUE) )
printME( b5effF <- efficiencies( b5, minusU = FALSE, margEff = TRUE ) )
printME( b5effDF <- efficiencies( b5, asInData = TRUE, minusU = FALSE, 
   margEff = TRUE ) )
b5m1 <- b5
b5m1$dataTable[ , "EDYRS" ] <- b5m1$dataTable[ , "EDYRS" ] + 1e-6
all.equal( attr( b5eff, "margEff" )[ , , 1 ], 
   ( efficiencies( b5m1 ) - b5eff )[ , ] / 1e-6, tol = 1e-4 )
all.equal( attr( b5effD, "margEff" )[ , 1 ], 
   c( efficiencies( b5m1, asInData = TRUE ) - b5effD ) / 1e-6, tol = 1e-4 )
all.equal( attr( b5effF, "margEff" )[ , , 1 ], 
   ( efficiencies( b5m1, minusU = FALSE ) - b5effF )[ ,  ] / 1e-6, tol = 1e-4 )
all.equal( attr( b5effDF, "margEff" )[ , 1 ],
   c( efficiencies( b5m1, asInData = TRUE, minusU = FALSE ) - b5effDF ) / 1e-6,
   tol = 1e-4 )
b5m2 <- b5
b5m2$dataTable[ , "BANRAT" ] <- b5m2$dataTable[ , "BANRAT" ] + 1e-6
all.equal( attr( b5eff, "margEff" )[ , , 2 ], 
   ( efficiencies( b5m2 ) - b5eff )[ , ] / 1e-6, tol = 1e-4 )
all.equal( attr( b5effD, "margEff" )[ , 2 ], 
   c( efficiencies( b5m2, asInData = TRUE ) - b5effD ) / 1e-6, tol = 1e-4 )
all.equal( attr( b5effF, "margEff" )[ , , 2 ], 
   ( efficiencies( b5m2, minusU = FALSE ) - b5effF )[ , ] / 1e-6, tol = 1e-4 )
all.equal( attr( b5effDF, "margEff" )[ , 2 ],
   c( efficiencies( b5m2, asInData = TRUE, minusU = FALSE ) - b5effDF ) / 1e-6,
   tol = 1e-4 )
all.equal( efficiencies( b5, margEff = TRUE ),
   efficiencies( sb5i, margEff = TRUE ) )
all.equal( efficiencies( b5, asInData = TRUE, margEff = TRUE, minusU = FALSE ), 
   efficiencies( sb5i, asInData = TRUE, margEff = TRUE, minusU = FALSE ) )
round( fitted( b5 ), 2 )
round( fitted( b5, asInData = TRUE ), 2 )
round( residuals( b5 ), 2 )
round( residuals( b5, asInData = TRUE ), 2 )
all.equal( fitted( b5, asInData = TRUE ) + residuals( b5, asInData = TRUE ),
   log( riceProdPhilPanel$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( b5 )

## panel data, efficiency effects frontier, zIntercept
sb6 <- sfa( lPROD ~ lAREA + lLABOR + lNPK | EDYRS + BANRAT,
   data = riceProdPhilPanel )
b6 <- frontier( data = riceProdPhilPanel,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   zNames = c( "EDYRS", "BANRAT" ), zIntercept = TRUE )
all.equal( sb6[-42], b6[-42], tol = 1e-4 )
sb6i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1 | EDYRS + BANRAT,
   data = riceProdPhilPanel )
all.equal( sb6i[ -c( 3, 7, 20, 42 ) ], sb6[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b6, digits = 1 )
coef( b6, which = "start" )
round( coef( b6, which = "ols" ), 2 )
round( coef( b6, which = "grid" ), 2 )
round( coef( b6, which = "grid", extraPar = TRUE ), 2 )
round( coef( b6 ), 2 )
round( coef( b6, extraPar = TRUE ), 2 )
round( coef( summary( b6 ), which = "ols" ), 2 )
round( coef( summary( b6 ) ), 2 )
round( coef( summary( b6, extraPar = TRUE ) ), 2 )
round( vcov( b6 ), 2 )
round( vcov( b6, extraPar = TRUE ), 2 )
print( logLik( b6, which = "ols" ), digits = 4 )
print( logLik( b6 ), digits = 4 )
nobs( b6 )
print( summary( b6 ), digits = 1 )
print( summary( b6, extraPar = TRUE ), digits = 1 )
all.equal( summary( b6, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb6i, extraPar = TRUE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( b6, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb6i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( b6 )
printME( efficiencies( b6, margEff = TRUE ) )
printME( efficiencies( b6, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( b6, margEff = TRUE, minusU = FALSE ),
   efficiencies( sb6i, margEff = TRUE, minusU = FALSE ) )
all.equal( efficiencies( b6, asInData = TRUE, margEff = TRUE ), 
   efficiencies( sb6i, asInData = TRUE, margEff = TRUE ) )
round( residuals( b6 ), 2 )
round( residuals( b6, asInData = TRUE ), 2 )
all.equal( fitted( b6, asInData = TRUE ) + residuals( b6, asInData = TRUE ),
   log( riceProdPhilPanel$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( b6 )

## panel data, error components frontier, truncNorm, timeEffect, starting values
sb7 <- sfa( lPROD ~ lAREA + lLABOR + lNPK, data = riceProdPhilPanel,
   truncNorm = TRUE, timeEffect = TRUE,
   startVal = c( -1, 0.3, 0.3, 0.3, 0.2, 0.5, -0.3, 0.1 ) )
b7 <- frontier( data = riceProdPhilPanel,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   truncNorm = TRUE, timeEffect = TRUE,
   startVal = c( -1, 0.3, 0.3, 0.3, 0.2, 0.5, -0.3, 0.1 ) )
all.equal( sb7[-42], b7[-42], tol = 1e-4 )
sb7i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1, 
   data = riceProdPhilPanel, truncNorm = TRUE, timeEffect = TRUE,
   startVal = c( -1, 0.3, 0.3, 0.3, 0.2, 0.5, -0.3, 0.1 ) )
all.equal( sb7i[ -c( 3, 7, 21, 42 ) ], sb7[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b7, digits = 1 )
coef( b7, which = "start" )
round( coef( b7, which = "ols" ), 2 )
coef( b7, which = "grid" )
coef( b7, which = "grid", extraPar = TRUE )
round( coef( b7 ), 2 )
round( coef( b7, extraPar = TRUE ), 2 )
round( coef( summary( b7 ), which = "ols" ), 2 )
round( coef( summary( b7 ) ), 2 )
round( coef( summary( b7, extraPar = TRUE ) ), 2 )
round( vcov( b7 ), 2 )
round( vcov( b7, extraPar = TRUE ), 2 )
print( logLik( b7, which = "ols" ), digits = 4 )
print( logLik( b7 ), digits = 4 )
nobs( b7 )
print( summary( b7 ), digits = 1 )
print( summary( b7, extraPar = TRUE ), digits = 1 )
all.equal( summary( b7, extraPar = TRUE )[ -c( 3, 7, 21, 42 ) ], 
   summary( sb7i, extraPar = TRUE )[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( b7, effMinusU = FALSE )[ -c( 3, 7, 21, 42 ) ], 
   summary( sb7i, effMinusU = FALSE )[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE )
lrtest( b7 )
round( efficiencies( b7 ), 2 )
round( efficiencies( b7, asInData = TRUE ), 2 )
all.equal( efficiencies( b7, minusU = FALSE ), 
   efficiencies( sb7i, minusU = FALSE ) )
all.equal( efficiencies( b7, asInData = TRUE ), 
   efficiencies( sb7i, asInData = TRUE ) )
all.equal( fitted( b7, asInData = TRUE ) + residuals( b7, asInData = TRUE ),
   log( riceProdPhilPanel$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( b7 )

## panel data, efficiency effects frontier, zIntercept, starting values
sb8 <- sfa( lPROD ~ lAREA + lLABOR + lNPK | EDYRS + BANRAT,
   data = riceProdPhilPanel,
   startVal = c( -1, 0.3, 0.3, 0.3, -0.3, -0.01, -0.4, 0.2, 0.8 ) )
b8 <- frontier( data = riceProdPhilPanel,
   yName = "lPROD", xNames = c( "lAREA", "lLABOR", "lNPK" ),
   zNames = c( "EDYRS", "BANRAT" ), zIntercept = TRUE,
   startVal = c( -1, 0.3, 0.3, 0.3, -0.3, -0.01, -0.4, 0.2, 0.8 ) )
all.equal( sb8[-42], b8[-42], tol = 1e-4 )
sb8i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1 | EDYRS + BANRAT, 
   data = riceProdPhilPanel, 
   startVal = c( -1, 0.3, 0.3, 0.3, -0.3, -0.01, -0.4, 0.2, 0.8 ) )
all.equal( sb8i[ -c( 3, 7, 21, 42 ) ], sb8[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b8, digits = 1 )
coef( b8, which = "start" )
round( coef( b8, which = "ols" ), 2 )
coef( b8, which = "grid" )
coef( b8, which = "grid", extraPar = TRUE )
round( coef( b8 ), 2 )
round( coef( b8, extraPar = TRUE ), 2 )
round( coef( summary( b8 ), which = "ols" ), 2 )
round( coef( summary( b8 ) ), 2 )
round( coef( summary( b8, extraPar = TRUE ) ), 2 )
round( vcov( b8 ), 2 )
round( vcov( b8, extraPar = TRUE ), 2 )
print( logLik( b8, which = "ols" ), digits = 4 )
print( logLik( b8 ), digits = 4 )
nobs( b8 )
print( summary( b8 ), digits = 1 )
print( summary( b8, extraPar = TRUE ), digits = 1 )
all.equal( summary( b8, extraPar = TRUE )[ -c( 3, 7, 21, 42 ) ], 
   summary( sb8i, extraPar = TRUE )[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE )
all.equal( summary( b8, effMinusU = FALSE )[ -c( 3, 7, 21, 42 ) ], 
   summary( sb8i, effMinusU = FALSE )[ -c( 3, 7, 21, 42 ) ], 
   check.attributes = FALSE )
lrtest( b8 )
printME( efficiencies( b8, margEff = TRUE ) )
printME( efficiencies( b8, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( b8 ), efficiencies( sb8i ) )
all.equal( efficiencies( b8, asInData = TRUE, minusU = FALSE ), 
   efficiencies( sb8i, asInData = TRUE, minusU = FALSE ) )
all.equal( fitted( b8, asInData = TRUE ) + residuals( b8, asInData = TRUE ),
   log( riceProdPhilPanel$PROD ), check.attributes = FALSE, tol = 1e-4 )
printAll( b8 )

## panel data, efficiency effects frontier: no Z vars
b9 <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ) | - 1,
   data = riceProdPhilPanel )
sb9i <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1 | - 1, 
   data = riceProdPhilPanel )
all.equal( sb9i[ -c( 3, 7, 20, 42 ) ], b9[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
round( coef( summary( b9, extraPar = TRUE ) ), 2 )
print( summary( b9 ), digits = 1 )
all.equal( summary( b9 )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb9i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b9, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sb9i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b9 )
lrtest( b9 )
round( efficiencies( b9 ), 2 )
round( efficiencies( b9, asInData = TRUE ), 2 )
all.equal( efficiencies( b9 ), efficiencies( sb9i ) )
all.equal( efficiencies( b9, asInData = TRUE, minusU = FALSE ), 
   efficiencies( sb9i, asInData = TRUE, minusU = FALSE ) )
all.equal( fitted( b9, asInData = TRUE ) + residuals( b9, asInData = TRUE ),
   log( riceProdPhilPanel$PROD ), check.attributes = FALSE, tol = 1e-4 )


## Cost Frontier (with land as quasi-fixed input)
## panel rice data, error components cost frontier
sd1 <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP, data = riceProdPhilPanel,
   ineffDecrease = FALSE )
Sd1 <- sfa( log( cost ) ~ log( PROD ) + log( AREA ) + log( LABORP ) +
   log( NPKP ), data = riceProdPhilPanel, ineffDecrease = FALSE )
all.equal( Sd1[-42], sd1[-42], check.attributes = FALSE, tol = 1e-4 )
d1 <- frontier( "lCost", xNames = c( "lPROD", "lAREA", "lLABORP", "lNPKP" ),
   data = riceProdPhilPanel, ineffDecrease = FALSE )
all.equal( sd1[-42], d1[-42], tol = 1e-4 )
sd1i <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1, 
   data = riceProdPhilPanel, ineffDecrease = FALSE )
all.equal( sd1i[ -c( 3, 7, 20, 42 ) ], sd1[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( d1, digits = 1 )
coef( d1, which = "start" )
round( coef( d1, which = "ols" ), 2 )
round( coef( d1, which = "grid" ), 2 )
round( coef( d1 ), 2 )
round( coef( summary( d1 ), which = "ols" ), 2 )
round( coef( summary( d1 ) ), 2 )
round( coef( summary( d1, extraPar = TRUE ) ), 2 )
round( vcov( d1 ), 2 )
nobs( d1 )
print( summary( d1 ), digits = 1 )
print( summary( d1, effMinusU = FALSE ), digits = 1 )
all.equal( summary( d1 )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd1i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( d1, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd1i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( d1 )
round( efficiencies( d1 ), 2 )
round( efficiencies( d1, asInData = TRUE ), 2 )
round( efficiencies( d1, minusU = FALSE ), 2 )
round( efficiencies( d1, asInData = TRUE, minusU = FALSE ), 2 )
all.equal( efficiencies( d1, minusU = FALSE ), 
   efficiencies( sd1i, minusU = FALSE ) )
all.equal( efficiencies( d1, asInData = TRUE ), 
   efficiencies( sd1i, asInData = TRUE ) )
round( residuals( d1 ), 2 )
round( residuals( d1, asInData = TRUE ), 2 )
all.equal( fitted( d1, asInData = TRUE ) + residuals( d1, asInData = TRUE ),
   log( riceProdPhilPanel$cost ), check.attributes = FALSE, tol = 1e-4 )
printAll( d1 )

## panel rice data, error components cost frontier, truncNorm
sd2 <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP, data = riceProdPhilPanel,
   ineffDecrease = FALSE, truncNorm = TRUE )
d2 <- frontier( "lCost", xNames = c( "lPROD", "lAREA", "lLABORP", "lNPKP" ),
   data = riceProdPhilPanel, ineffDecrease = FALSE, truncNorm = TRUE )
all.equal( sd2[-42], d2[-42], tol = 1e-4 )
sd2i <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1, 
   data = riceProdPhilPanel, ineffDecrease = FALSE, truncNorm = TRUE )
all.equal( sd2i[ -c( 3, 7, 20, 42 ) ], sd2[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( d2, digits = 1 )
coef( d2, which = "start" )
round( coef( d2, which = "ols" ), 2 )
round( coef( d2, which = "grid" ), 2 )
round( coef( d2 ), 2 )
round( coef( summary( d2 ), which = "ols" ), 2 )
round( coef( summary( d2 ) ), 2 )
round( coef( summary( d2, extraPar = TRUE ) ), 2 )
round( vcov( d2 ), 2 )
nobs( d2 )
print( summary( d2, effMinusU = FALSE ), digits = 1 )
all.equal( summary( d2 )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd2i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( d2, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd2i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( d2 )
round( efficiencies( d2, minusU = FALSE ), 2 )
round( efficiencies( d2, asInData = TRUE, minusU = FALSE ), 2 )
all.equal( efficiencies( d2 ), efficiencies( sd2i ) )
all.equal( efficiencies( d2, asInData = TRUE, minusU = FALSE ), 
   efficiencies( sd2i, asInData = TRUE, minusU = FALSE ) )
round( fitted( d2 ), 2 )
round( fitted( d2, asInData = TRUE ), 2 )
round( residuals( d2 ), 2 )
round( residuals( d2, asInData = TRUE ), 2 )
all.equal( fitted( d2, asInData = TRUE ) + residuals( d2, asInData = TRUE ),
   log( riceProdPhilPanel$cost ), check.attributes = FALSE, tol = 1e-4 )
printAll( d2 )

## panel rice data, error components cost frontier, timeEffect
sd3 <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP, data = riceProdPhilPanel,
   ineffDecrease = FALSE, timeEffect = TRUE )
d3 <- frontier( "lCost", xNames = c( "lPROD", "lAREA", "lLABORP", "lNPKP" ),
   data = riceProdPhilPanel, ineffDecrease = FALSE, timeEffect = TRUE )
all.equal( sd3[-42], d3[-42], tol = 1e-4 )
sd3i <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1, 
   data = riceProdPhilPanel, ineffDecrease = FALSE, timeEffect = TRUE )
all.equal( sd3i[ -c( 3, 7, 20, 42 ) ], sd3[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( d3, digits = 1 )
coef( d3, which = "start" )
round( coef( d3, which = "ols" ), 2 )
round( coef( d3, which = "grid" ), 2 )
round( coef( d3 ), 2 )
round( coef( summary( d3 ), which = "ols" ), 2 )
round( coef( summary( d3 ) ), 2 )
round( coef( summary( d3, extraPar = TRUE ) ), 2 )
round( vcov( d3 ), 2 )
nobs( d3 )
print( summary( d3, effMinusU = FALSE ), digits = 1 )
all.equal( summary( d3 )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd3i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( d3, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd3i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( d3 )
round( efficiencies( d3, minusU = FALSE ), 2 )
round( efficiencies( d3, asInData = TRUE, minusU = FALSE ), 2 )
all.equal( efficiencies( d3, minusU = FALSE ), 
   efficiencies( sd3i, minusU = FALSE ) )
all.equal( efficiencies( d3, asInData = TRUE ), 
   efficiencies( sd3i, asInData = TRUE ) )
round( residuals( d3 ), 2 )
round( residuals( d3, asInData = TRUE ), 2 )
all.equal( fitted( d3, asInData = TRUE ) + residuals( d3, asInData = TRUE ),
   log( riceProdPhilPanel$cost ), check.attributes = FALSE, tol = 1e-4 )
printAll( d3 )

## panel rice data, error components cost frontier, truncNorm, timeEffect
sd4 <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP, data = riceProdPhilPanel,
   ineffDecrease = FALSE, truncNorm = TRUE, timeEffect = TRUE )
d4 <- frontier( "lCost", xNames = c( "lPROD", "lAREA", "lLABORP", "lNPKP" ),
   data = riceProdPhilPanel, ineffDecrease = FALSE, truncNorm = TRUE,
   timeEffect = TRUE )
all.equal( sd4[-42], d4[-42], tol = 1e-4 )
sd4i <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1, 
   data = riceProdPhilPanel, ineffDecrease = FALSE, 
   truncNorm = TRUE, timeEffect = TRUE )
all.equal( sd4i[ -c( 3, 7, 20, 42 ) ], sd4[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( d4, digits = 1 )
coef( d4, which = "start" )
round( coef( d4, which = "ols" ), 2 )
round( coef( d4, which = "grid" ), 2 )
round( coef( d4 ), 2 )
round( coef( summary( d4 ), which = "ols" ), 2 )
round( coef( summary( d4 ) ), 2 )
round( coef( summary( d4, extraPar = TRUE ) ), 2 )
round( vcov( d4 ), 2 )
nobs( d4 )
print( summary( d4, effMinusU = FALSE ), digits = 1 )
all.equal( summary( d4 )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd4i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( d4, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd4i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( d4 )
round( efficiencies( d4, minusU = FALSE ), 2 )
round( efficiencies( d4, asInData = TRUE, minusU = FALSE ), 2 )
all.equal( efficiencies( d4 ), efficiencies( sd4i ) )
all.equal( efficiencies( d4, asInData = TRUE, minusU = FALSE ), 
   efficiencies( sd4i, asInData = TRUE, minusU = FALSE ) )
round( residuals( d4 ), 2 )
round( residuals( d4, asInData = TRUE ), 2 )
all.equal( fitted( d4, asInData = TRUE ) + residuals( d4, asInData = TRUE ),
   log( riceProdPhilPanel$cost ), check.attributes = FALSE, tol = 1e-4 )
printAll( d4 )

## panel rice data, efficiency effects cost frontier
sd5 <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP | EDYRS + BANRAT - 1,
   data = riceProdPhilPanel, ineffDecrease = FALSE )
Sd5 <- sfa( log( cost ) ~ log( PROD ) + log( AREA ) + log( LABORP ) +
   log( NPKP ) | EDYRS + BANRAT - 1, data = riceProdPhilPanel,
   ineffDecrease = FALSE )
all.equal( Sd5[-42], sd5[-42], check.attributes = FALSE, tol = 1e-4 )
d5 <- frontier( "lCost", xNames = c( "lPROD", "lAREA", "lLABORP", "lNPKP" ),
   zNames = c( "EDYRS", "BANRAT" ), data = riceProdPhilPanel,
   ineffDecrease = FALSE )
all.equal( sd5[-42], d5[-42], tol = 1e-4 )
sd5i <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1 | 
      EDYRS + BANRAT - 1, data = riceProdPhilPanel, ineffDecrease = FALSE )
all.equal( sd5i[ -c( 3, 7, 20, 42 ) ], sd5[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( d5, digits = 1 )
coef( d5, which = "start" )
round( coef( d5, which = "ols" ), 2 )
round( coef( d5, which = "grid" ), 2 )
round( coef( d5 ), 2 )
round( coef( summary( d5 ), which = "ols" ), 2 )
round( coef( summary( d5 ) ), 2 )
round( coef( summary( d5, extraPar = TRUE ) ), 2 )
round( vcov( d5 ), 2 )
nobs( d5 )
print( summary( d5 ), digits = 1 )
print( summary( d5, effMinusU = FALSE ), digits = 1 )
all.equal( summary( d5 )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd5i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( d5, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd5i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( d5 )
printME( d5eff <- efficiencies( d5, margEff = TRUE ) )
printME( d5effD <- efficiencies( d5, asInData = TRUE, margEff = TRUE ) )
printME( d5effF <- efficiencies( d5, minusU = FALSE, margEff = TRUE ) )
printME( d5effDF <- efficiencies( d5, asInData = TRUE, minusU = FALSE, 
   margEff = TRUE ) )
d5m1 <- d5
d5m1$dataTable[ , "EDYRS" ] <- d5m1$dataTable[ , "EDYRS" ] + 1e-6
all.equal( attr( d5eff, "margEff" )[ , , 1 ], 
   ( efficiencies( d5m1 ) - d5eff )[ , ] / 1e-6, tol = 1e-4 )
all.equal( attr( d5effD, "margEff" )[ , 1 ], 
   c( efficiencies( d5m1, asInData = TRUE ) - d5effD ) / 1e-6, tol = 1e-4 )
all.equal( attr( d5effF, "margEff" )[ , , 1 ], 
   ( efficiencies( d5m1, minusU = FALSE ) - d5effF )[ ,  ] / 1e-6, tol = 1e-4 )
all.equal( attr( d5effDF, "margEff" )[ , 1 ],
   c( efficiencies( d5m1, asInData = TRUE, minusU = FALSE ) - d5effDF ) / 1e-6,
   tol = 1e-4 )
d5m2 <- d5
d5m2$dataTable[ , "BANRAT" ] <- d5m2$dataTable[ , "BANRAT" ] + 1e-6
all.equal( attr( d5eff, "margEff" )[ , , 2 ], 
   ( efficiencies( d5m2 ) - d5eff )[ , ] / 1e-6, tol = 1e-4 )
all.equal( attr( d5effD, "margEff" )[ , 2 ], 
   c( efficiencies( d5m2, asInData = TRUE ) - d5effD ) / 1e-6, tol = 1e-4 )
all.equal( attr( d5effF, "margEff" )[ , , 2 ], 
   ( efficiencies( d5m2, minusU = FALSE ) - d5effF )[ , ] / 1e-6, tol = 1e-4 )
all.equal( attr( d5effDF, "margEff" )[ , 2 ],
   c( efficiencies( d5m2, asInData = TRUE, minusU = FALSE ) - d5effDF ) / 1e-6,
   tol = 1e-4 )
all.equal( efficiencies( d5, margEff = TRUE, minusU = FALSE ), 
   efficiencies( sd5i, margEff = TRUE, minusU = FALSE ) )
all.equal( efficiencies( d5, asInData = TRUE, margEff = TRUE ), 
   efficiencies( sd5i, asInData = TRUE, margEff = TRUE ) )
round( fitted( d5 ), 2 )
round( fitted( d5, asInData = TRUE ), 2 )
round( residuals( d5 ), 2 )
round( residuals( d5, asInData = TRUE ), 2 )
all.equal( fitted( d5, asInData = TRUE ) + residuals( d5, asInData = TRUE ),
   log( riceProdPhilPanel$cost ), check.attributes = FALSE, tol = 1e-4 )
printAll( d5 )

## panel rice data, efficiency effects cost frontier, zIntercept
sd6 <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP | EDYRS + BANRAT,
   data = riceProdPhilPanel, ineffDecrease = FALSE )
d6 <- frontier( "lCost", xNames = c( "lPROD", "lAREA", "lLABORP", "lNPKP" ),
   zNames = c( "EDYRS", "BANRAT" ), data = riceProdPhilPanel,
   ineffDecrease = FALSE, zIntercept = TRUE )
all.equal( sd6[-42], d6[-42], tol = 1e-4 )
sd6i <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1 | 
      EDYRS + BANRAT, data = riceProdPhilPanel, ineffDecrease = FALSE )
all.equal( sd6i[ -c( 3, 7, 20, 42 ) ], sd6[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( d6, digits = 1 )
coef( d6, which = "start" )
round( coef( d6, which = "ols" ), 2 )
round( coef( d6, which = "grid" ), 2 )
round( coef( d6 ), 2 )
round( coef( summary( d6 ), which = "ols" ), 2 )
round( coef( summary( d6 ) ), 2 )
round( coef( summary( d6, extraPar = TRUE ) ), 2 )
round( vcov( d6 ), 2 )
nobs( d6 )
print( summary( d6, effMinusU = FALSE ), digits = 1 )
all.equal( summary( d6 )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd6i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( d6, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd6i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( d6 )
printME( efficiencies( d6, minusU = FALSE, margEff = TRUE ) )
printME( efficiencies( d6, asInData = TRUE, minusU = FALSE, margEff = TRUE ) )
all.equal( efficiencies( d6, margEff = TRUE ), 
   efficiencies( sd6i, margEff = TRUE ) )
all.equal( efficiencies( d6, asInData = TRUE, margEff = TRUE, minusU = FALSE ), 
   efficiencies( sd6i, asInData = TRUE, margEff = TRUE, minusU = FALSE ) )
round( residuals( d6 ), 2 )
round( residuals( d6, asInData = TRUE ), 2 )
all.equal( fitted( d6, asInData = TRUE ) + residuals( d6, asInData = TRUE ),
   log( riceProdPhilPanel$cost ), check.attributes = FALSE, tol = 1e-4 )
printAll( d6 )

## panel rice data, efficiency effects cost frontier: no Z vars
d9 <- sfa( log( cost ) ~ log( PROD ) + log( AREA ) + log( LABORP ) +
   log( NPKP ) | - 1, data = riceProdPhilPanel, ineffDecrease = FALSE )
sd9i <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1 | - 1, 
   data = riceProdPhilPanel, ineffDecrease = FALSE )
all.equal( sd9i[ -c( 3, 7, 20, 42 ) ], d9[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
round( coef( summary( d9, extraPar = TRUE ) ), 2 )
print( summary( d9, effMinusU = FALSE ), digits = 1 )
all.equal( summary( d9 )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd9i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( d9, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( sd9i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( d9 )
lrtest( d9 )
round( efficiencies( d9, minusU = FALSE ), 2 )
round( efficiencies( d9, asInData = TRUE ), 2 )
all.equal( efficiencies( d9, minusU = FALSE ), 
   efficiencies( sd9i, minusU = FALSE ) )
all.equal( efficiencies( d9, asInData = TRUE ), 
   efficiencies( sd9i, asInData = TRUE ) )
all.equal( fitted( d9, asInData = TRUE ) + residuals( d9, asInData = TRUE ),
   log( riceProdPhilPanel$cost ), check.attributes = FALSE, tol = 1e-4 )


## unbalanced panel data
set.seed( 321 )
riceProdPhilPanelUnb <- riceProdPhilPanel
riceProdPhilPanelUnb[ 3, c( "PROD", "lPROD" ) ] <- NA
riceProdPhilPanelUnb[ 5, c( "AREA", "lAREA" ) ] <- NA
riceProdPhilPanelUnb[ 111, c( "LABOR", "lLABOR", "LABORP", "lLABORP" ) ] <- NA
riceProdPhilPanelUnb[ 222, c( "NPK", "lNPK", "NPKP", "lNPKP" ) ] <- NA

## unbalanced panel data, error components frontier
b1u <- sfa( lPROD ~ lAREA + lLABOR + lNPK, data = riceProdPhilPanelUnb )
b1ui <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1, 
   data = riceProdPhilPanelUnb )
all.equal( b1ui[ -c( 3, 7, 20, 42 ) ], b1u[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b1u, digits = 1 )
print( summary( b1u ), digits = 1 )
all.equal( summary( b1u )[ -c( 3, 7, 20, 42 ) ], 
   summary( b1ui )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b1u, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b1ui, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b1u )
lrtest( b1u )
round( efficiencies( b1u ), 2 )
round( efficiencies( b1u, asInData = TRUE ), 2 )
all.equal( efficiencies( b1u ), efficiencies( b1ui ) )
all.equal( efficiencies( b1u, asInData = TRUE, minusU = FALSE ), 
   efficiencies( b1ui, asInData = TRUE, minusU = FALSE ) )
round( fitted( b1u ), 2 )
round( fitted( b1u, asInData = TRUE ), 2 )
round( residuals( b1u ), 2 )
round( residuals( b1u, asInData = TRUE ), 2 )
all.equal( fitted( b1u, asInData = TRUE ) + residuals( b1u, asInData = TRUE ),
   log( riceProdPhilPanelUnb$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( riceProdPhilPanelUnb$PROD ) - fitted( b1u, asInData = TRUE ),
   residuals( b1u, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b1u )

## unbalanced panel data, error components frontier, truncNorm
b2u <- sfa( lPROD ~ lAREA + lLABOR + lNPK, data = riceProdPhilPanelUnb,
   truncNorm = TRUE )
b2ui <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1, 
   data = riceProdPhilPanelUnb, truncNorm = TRUE )
all.equal( b2ui[ -c( 3, 7, 20, 42 ) ], b2u[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b2u, digits = 1 )
print( summary( b2u ), digits = 1 )
all.equal( summary( b2u )[ -c( 3, 7, 20, 42 ) ], 
   summary( b2ui )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b2u, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b2ui, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b2u )
lrtest( b2u )
round( efficiencies( b2u ), 2 )
round( efficiencies( b2u, asInData = TRUE ), 2 )
all.equal( efficiencies( b2u ), efficiencies( b2ui ) )
all.equal( efficiencies( b2u, asInData = TRUE, minusU = FALSE ), 
   efficiencies( b2ui, asInData = TRUE, minusU = FALSE ) )
round( residuals( b2u ), 2 )
round( residuals( b2u, asInData = TRUE ), 2 )
all.equal( fitted( b2u, asInData = TRUE ) + residuals( b2u, asInData = TRUE ),
   log( riceProdPhilPanelUnb$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( riceProdPhilPanelUnb$PROD ) - fitted( b2u, asInData = TRUE ),
   residuals( b2u, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b2u )

## unbalanced panel data, error components frontier, timeEffect
b3u <- sfa( lPROD ~ lAREA + lLABOR + lNPK, data = riceProdPhilPanelUnb,
   timeEffect = TRUE )
b3ui <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1, 
   data = riceProdPhilPanelUnb, timeEffect = TRUE )
all.equal( b3ui[ -c( 3, 7, 20, 42 ) ], b3u[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b3u, digits = 1 )
print( summary( b3u ), digits = 1 )
all.equal( summary( b3u )[ -c( 3, 7, 20, 42 ) ], 
   summary( b3ui )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b3u, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b3ui, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b3u )
lrtest( b3u )
round( efficiencies( b3u ), 2 )
round( efficiencies( b3u, asInData = TRUE ), 2 )
all.equal( efficiencies( b3u, minusU = FALSE ), 
   efficiencies( b3ui, minusU = FALSE ) )
all.equal( efficiencies( b3u, asInData = TRUE ), 
   efficiencies( b3ui, asInData = TRUE ) )
round( residuals( b3u ), 2 )
round( residuals( b3u, asInData = TRUE ), 2 )
all.equal( fitted( b3u, asInData = TRUE ) + residuals( b3u, asInData = TRUE ),
   log( riceProdPhilPanelUnb$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( riceProdPhilPanelUnb$PROD ) - fitted( b3u, asInData = TRUE ),
   residuals( b3u, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b3u )

## unbalanced panel data, error components frontier, truncNorm, timeEffect
b4u <- sfa( lPROD ~ lAREA + lLABOR + lNPK, data = riceProdPhilPanelUnb,
   truncNorm = TRUE, timeEffect = TRUE )
b4ui <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1, 
   data = riceProdPhilPanelUnb, truncNorm = TRUE, timeEffect = TRUE )
all.equal( b4ui[ -c( 3, 7, 20, 42 ) ], b4u[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b4u, digits = 1 )
print( summary( b4u ), digits = 1 )
print( summary( b4u, effMinusU = FALSE ), digits = 1 )
all.equal( summary( b4u )[ -c( 3, 7, 20, 42 ) ], 
   summary( b4ui )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b4u, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b4ui, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b4u )
lrtest( b4u )
round( efficiencies( b4u ), 2 )
round( efficiencies( b4u, asInData = TRUE ), 2 )
round( efficiencies( b4u, minusU = FALSE ), 2 )
round( efficiencies( b4u, asInData = TRUE, minusU = FALSE ), 2 )
all.equal( efficiencies( b4u ), efficiencies( b4ui ) )
all.equal( efficiencies( b4u, asInData = TRUE, minusU = FALSE ), 
   efficiencies( b4ui, asInData = TRUE, minusU = FALSE ) )
round( residuals( b4u ), 2 )
round( residuals( b4u, asInData = TRUE ), 2 )
all.equal( fitted( b4u, asInData = TRUE ) + residuals( b4u, asInData = TRUE ),
   log( riceProdPhilPanelUnb$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( riceProdPhilPanelUnb$PROD ) - fitted( b4u, asInData = TRUE ),
   residuals( b4u, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b4u )

## unbalanced panel data, efficiency effects frontier
b5u <- sfa( lPROD ~ lAREA + lLABOR + lNPK | EDYRS + BANRAT - 1,
   data = riceProdPhilPanelUnb )
b5ui <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1 | EDYRS + BANRAT - 1, 
   data = riceProdPhilPanelUnb )
all.equal( b5ui[ -c( 3, 7, 20, 42 ) ], b5u[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b5u, digits = 1 )
print( summary( b5u ), digits = 1 )
all.equal( summary( b5u )[ -c( 3, 7, 20, 42 ) ], 
   summary( b5ui )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b5u, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b5ui, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b5u )
lrtest( b5u )
printME( efficiencies( b5u, margEff = TRUE ) )
printME( efficiencies( b5u, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( b5u, margEff = TRUE, minusU = FALSE ), 
   efficiencies( b5ui, margEff = TRUE, minusU = FALSE ) )
all.equal( efficiencies( b5u, asInData = TRUE, margEff = TRUE ), 
   efficiencies( b5ui, asInData = TRUE, margEff = TRUE ) )
round( residuals( b5u ), 2 )
round( residuals( b5u, asInData = TRUE ), 2 )
all.equal( fitted( b5u, asInData = TRUE ) + residuals( b5u, asInData = TRUE ),
   log( riceProdPhilPanelUnb$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( riceProdPhilPanelUnb$PROD ) - fitted( b5u, asInData = TRUE ),
   residuals( b5u, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b5u )

## unbalanced panel data, efficiency effects frontier, zIntercept
b6u <- sfa( lPROD ~ lAREA + lLABOR + lNPK | EDYRS + BANRAT,
   data = riceProdPhilPanelUnb )
b6ui <- sfa( lPROD ~ ones + lAREA + lLABOR + lNPK - 1 | EDYRS + BANRAT, 
   data = riceProdPhilPanelUnb )
all.equal( b6ui[ -c( 3, 7, 20, 42 ) ], b6u[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b6u, digits = 1 )
print( summary( b6u ), digits = 1 )
print( summary( b6u, effMinusU = FALSE ), digits = 1 )
all.equal( summary( b6u )[ -c( 3, 7, 20, 42 ) ], 
   summary( b6ui )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b6u, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b6ui, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b6u )
lrtest( b6u )
printME( b6ueff <- efficiencies( b6u, margEff = TRUE ) )
printME( b6ueffD <- efficiencies( b6u, asInData = TRUE, margEff = TRUE ) )
printME( b6ueffF <- efficiencies( b6u, minusU = FALSE, margEff = TRUE ) )
printME( b6ueffDF <- efficiencies( b6u, asInData = TRUE, minusU = FALSE, 
   margEff = TRUE ) )
b6um1 <- b6u
b6um1$dataTable[ , "EDYRS" ] <- b6um1$dataTable[ , "EDYRS" ] + 1e-6
all.equal( attr( b6ueff, "margEff" )[ , , 1 ], 
   ( efficiencies( b6um1 ) - b6ueff )[ , ] / 1e-6, tol = 1e-4 )
all.equal( attr( b6ueffD, "margEff" )[ , 1 ], 
   c( efficiencies( b6um1, asInData = TRUE ) - b6ueffD ) / 1e-6, tol = 1e-4 )
all.equal( attr( b6ueffF, "margEff" )[ , , 1 ], 
   ( efficiencies( b6um1, minusU = FALSE ) - b6ueffF )[ ,  ] / 1e-6, tol = 1e-4 )
all.equal( attr( b6ueffDF, "margEff" )[ , 1 ],
   c( efficiencies( b6um1, asInData = TRUE, minusU = FALSE ) - b6ueffDF ) / 1e-6,
   tol = 1e-4 )
b6um2 <- b6u
b6um2$dataTable[ , "BANRAT" ] <- b6um2$dataTable[ , "BANRAT" ] + 1e-6
all.equal( attr( b6ueff, "margEff" )[ , , 2 ], 
   ( efficiencies( b6um2 ) - b6ueff )[ , ] / 1e-6, tol = 1e-4 )
all.equal( attr( b6ueffD, "margEff" )[ , 2 ], 
   c( efficiencies( b6um2, asInData = TRUE ) - b6ueffD ) / 1e-6, tol = 1e-4 )
all.equal( attr( b6ueffF, "margEff" )[ , , 2 ], 
   ( efficiencies( b6um2, minusU = FALSE ) - b6ueffF )[ , ] / 1e-6, tol = 1e-4 )
all.equal( attr( b6ueffDF, "margEff" )[ , 2 ],
   c( efficiencies( b6um2, asInData = TRUE, minusU = FALSE ) - b6ueffDF ) / 1e-6,
   tol = 1e-4 )
all.equal( efficiencies( b6u, margEff = TRUE ), 
   efficiencies( b6ui, margEff = TRUE ) )
all.equal( efficiencies( b6u, asInData = TRUE, margEff = TRUE, minusU = FALSE ), 
   efficiencies( b6ui, asInData = TRUE, margEff = TRUE, minusU = FALSE ) )
round( fitted( b6u ), 2 )
round( fitted( b6u, asInData = TRUE ), 2 )
round( residuals( b6u ), 2 )
round( residuals( b6u, asInData = TRUE ), 2 )
all.equal( fitted( b6u, asInData = TRUE ) + residuals( b6u, asInData = TRUE ),
   log( riceProdPhilPanelUnb$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( riceProdPhilPanelUnb$PROD ) - fitted( b6u, asInData = TRUE ),
   residuals( b6u, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b6u )

## unbalanced panel rice data, error components cost frontier
d1u <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP, data = riceProdPhilPanelUnb,
   ineffDecrease = FALSE )
d1ui <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1, 
   data = riceProdPhilPanelUnb, ineffDecrease = FALSE )
all.equal( d1ui[ -c( 3, 7, 20, 42 ) ], d1u[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( d1u, digits = 1 )
print( summary( d1u ), digits = 1 )
print( summary( d1u, effMinusU = FALSE ), digits = 1 )
all.equal( summary( d1u )[ -c( 3, 7, 20, 42 ) ], 
   summary( d1ui )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( d1u, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( d1ui, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( d1u )
lrtest( d1u )
round( efficiencies( d1u ), 2 )
round( efficiencies( d1u, asInData = TRUE ), 2 )
round( efficiencies( d1u, minusU = FALSE ), 2 )
round( efficiencies( d1u, asInData = TRUE, minusU = FALSE ), 2 )
all.equal( efficiencies( d1u ), efficiencies( d1ui ) )
all.equal( efficiencies( d1u, asInData = TRUE, minusU = FALSE ), 
   efficiencies( d1ui, asInData = TRUE, minusU = FALSE ) )
round( residuals( d1u ), 2 )
round( residuals( d1u, asInData = TRUE ), 2 )
all.equal( fitted( d1u, asInData = TRUE ) + residuals( d1u, asInData = TRUE ),
   log( riceProdPhilPanelUnb$cost ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( riceProdPhilPanelUnb$cost ) - fitted( d1u, asInData = TRUE ),
   residuals( d1u, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( d1u )

## unbalanced panel rice data, error components cost frontier, truncNorm
d2u <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP, data = riceProdPhilPanelUnb,
   ineffDecrease = FALSE, truncNorm = TRUE )
d2ui <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1, 
   data = riceProdPhilPanelUnb, ineffDecrease = FALSE, truncNorm = TRUE )
all.equal( d2ui[ -c( 3, 7, 20, 42 ) ], d2u[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( d2u, digits = 1 )
print( summary( d2u, effMinusU = FALSE ), digits = 1 )
all.equal( summary( d2u )[ -c( 3, 7, 20, 42 ) ], 
   summary( d2ui )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( d2u, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( d2ui, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( d2u )
lrtest( d2u )
round( efficiencies( d2u, minusU = FALSE ), 2 )
round( efficiencies( d2u, asInData = TRUE, minusU = FALSE ), 2 )
all.equal( efficiencies( d2u, minusU = FALSE ), 
   efficiencies( d2ui, minusU = FALSE ) )
all.equal( efficiencies( d2u, asInData = TRUE ), 
   efficiencies( d2ui, asInData = TRUE ) )
round( residuals( d2u ), 2 )
round( residuals( d2u, asInData = TRUE ), 2 )
all.equal( fitted( d2u, asInData = TRUE ) + residuals( d2u, asInData = TRUE ),
   log( riceProdPhilPanelUnb$cost ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( riceProdPhilPanelUnb$cost ) - fitted( d2u, asInData = TRUE ),
   residuals( d2u, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( d2u )

## unbalanced panel rice data, error components cost frontier, timeEffect
d3u <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP, data = riceProdPhilPanelUnb,
   ineffDecrease = FALSE, timeEffect = TRUE )
d3ui <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1, 
   data = riceProdPhilPanelUnb, ineffDecrease = FALSE, timeEffect = TRUE )
all.equal( d3ui[ -c( 3, 7, 20, 42 ) ], d3u[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( d3u, digits = 1 )
print( summary( d3u, effMinusU = FALSE ), digits = 1 )
all.equal( summary( d3u )[ -c( 3, 7, 20, 42 ) ], 
   summary( d3ui )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( d3u, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( d3ui, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( d3u )
lrtest( d3u )
round( efficiencies( d3u, minusU = FALSE ), 2 )
round( efficiencies( d3u, asInData = TRUE, minusU = FALSE ), 2 )
all.equal( efficiencies( d3u ), efficiencies( d3ui ) )
all.equal( efficiencies( d3u, asInData = TRUE, minusU = FALSE ), 
   efficiencies( d3ui, asInData = TRUE, minusU = FALSE ) )
round( fitted( d3u ), 2 )
round( fitted( d3u, asInData = TRUE ), 2 )
round( residuals( d3u ), 2 )
round( residuals( d3u, asInData = TRUE ), 2 )
all.equal( fitted( d3u, asInData = TRUE ) + residuals( d3u, asInData = TRUE ),
   log( riceProdPhilPanelUnb$cost ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( riceProdPhilPanelUnb$cost ) - fitted( d3u, asInData = TRUE ),
   residuals( d3u, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( d3u )

## unbalanced panel rice data, error components cost frontier, truncNorm, timeEffect
d4u <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP, data = riceProdPhilPanelUnb,
   ineffDecrease = FALSE, truncNorm = TRUE, timeEffect = TRUE )
d4ui <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1, 
   data = riceProdPhilPanelUnb, ineffDecrease = FALSE, 
   truncNorm = TRUE, timeEffect = TRUE )
all.equal( d4ui[ -c( 3, 7, 20, 42 ) ], d4u[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( d4u, digits = 1 )
print( summary( d4u, effMinusU = FALSE ), digits = 1 )
all.equal( summary( d4u )[ -c( 3, 7, 20, 42 ) ], 
   summary( d4ui )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( d4u, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( d4ui, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( d4u )
lrtest( d4u )
round( efficiencies( d4u, minusU = FALSE ), 2 )
round( efficiencies( d4u, asInData = TRUE, minusU = FALSE ), 2 )
all.equal( efficiencies( d4u, minusU = FALSE ), 
   efficiencies( d4ui, minusU = FALSE ) )
all.equal( efficiencies( d4u, asInData = TRUE ), 
   efficiencies( d4ui, asInData = TRUE ) )
round( residuals( d4u ), 2 )
round( residuals( d4u, asInData = TRUE ), 2 )
all.equal( fitted( d4u, asInData = TRUE ) + residuals( d4u, asInData = TRUE ),
   log( riceProdPhilPanelUnb$cost ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( riceProdPhilPanelUnb$cost ) - fitted( d4u, asInData = TRUE ),
   residuals( d4u, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( d4u )

## unbalanced panel rice data, efficiency effects cost frontier
d5u <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP | EDYRS + BANRAT - 1,
   data = riceProdPhilPanelUnb, ineffDecrease = FALSE )
d5ui <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1 | 
      EDYRS + BANRAT - 1, data = riceProdPhilPanelUnb, ineffDecrease = FALSE )
all.equal( d5ui[ -c( 3, 7, 20, 42 ) ], d5u[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( d5u, minusU = FALSE )
print( summary( d5u ), digits = 1 )
print( summary( d5u, effMinusU = FALSE ), digits = 1 )
all.equal( summary( d5u )[ -c( 3, 7, 20, 42 ) ], 
   summary( d5ui )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( d5u, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( d5ui, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( d5u )
lrtest( d5u )
printME( efficiencies( d5u, margEff = TRUE ) )
printME( efficiencies( d5u, asInData = TRUE, margEff = TRUE ) )
printME( efficiencies( d5u, minusU = FALSE, margEff = TRUE ) )
printME( efficiencies( d5u, asInData = TRUE, minusU = FALSE, margEff = TRUE ) )
all.equal( efficiencies( d5u, margEff = TRUE ), 
   efficiencies( d5ui, margEff = TRUE ) )
all.equal( efficiencies( d5u, asInData = TRUE, margEff = TRUE, minusU = FALSE ), 
   efficiencies( d5ui, asInData = TRUE, margEff = TRUE, minusU = FALSE ) )
round( fitted( d5u ), 2 )
round( fitted( d5u, asInData = TRUE ), 2 )
round( residuals( d5u ), 2 )
round( residuals( d5u, asInData = TRUE ), 2 )
all.equal( fitted( d5u, asInData = TRUE ) + residuals( d5u, asInData = TRUE ),
   log( riceProdPhilPanelUnb$cost ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( riceProdPhilPanelUnb$cost ) - fitted( d5u, asInData = TRUE ),
   residuals( d5u, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( d5u )

## unbalanced panel rice data, efficiency effects cost frontier, zIntercept
d6u <- sfa( lCost ~ lPROD + lAREA + lLABORP + lNPKP | EDYRS + BANRAT,
   data = riceProdPhilPanelUnb, ineffDecrease = FALSE )
d6ui <- sfa( lCost ~ ones + lPROD + lAREA + lLABORP + lNPKP - 1 | 
      EDYRS + BANRAT, data = riceProdPhilPanelUnb, ineffDecrease = FALSE )
all.equal( d6ui[ -c( 3, 7, 20, 42 ) ], d6u[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( d6u, digits = 1 )
print( summary( d6u, effMinusU = FALSE ), digits = 1 )
all.equal( summary( d6u )[ -c( 3, 7, 20, 42 ) ], 
   summary( d6ui )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( d6u, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( d6ui, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( d6u )
lrtest( d6u )
printME( efficiencies( d6u, minusU = FALSE, margEff = TRUE ) )
printME( efficiencies( d6u, asInData = TRUE, minusU = FALSE, margEff = TRUE ) )
all.equal( efficiencies( d6u, margEff = TRUE, minusU = FALSE ), 
   efficiencies( d6ui, margEff = TRUE, minusU = FALSE ) )
all.equal( efficiencies( d6u, asInData = TRUE, margEff = TRUE ), 
   efficiencies( d6ui, asInData = TRUE, margEff = TRUE ) )
round( residuals( d6u ), 2 )
round( residuals( d6u, asInData = TRUE ), 2 )
all.equal( fitted( d6u, asInData = TRUE ) + residuals( d6u, asInData = TRUE ),
   log( riceProdPhilPanelUnb$cost ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( riceProdPhilPanelUnb$cost ) - fitted( d6u, asInData = TRUE ),
   residuals( d6u, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( d6u )


## unbalanced panel data with firms that have NAs in all time periods
naPanelData <- riceProdPhilPanelUnb
naPanelData[ naPanelData$farm == "F_21", "PROD" ] <- NA
naPanelData[ naPanelData$farm == "F_23", "AREA" ] <- NA
naPanelData[ naPanelData$farm == "F_26", "LABOR" ] <- NA
naPanelData[ naPanelData$farm == "F_30", "NPK" ] <- NA
naPanelData[ naPanelData$farm == "F_35", "EDYRS" ] <- NA

## panel data with NA firms, error components frontier
b1n <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ),
   data = naPanelData )
b1ni <- sfa( log( PROD ) ~ ones + log( AREA ) + log( LABOR ) + log( NPK ) - 1,
   data = naPanelData )
all.equal( b1ni[ -c( 3, 7, 20, 42 ) ], b1n[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b1n, digits = 1 )
print( summary( b1n ), digits = 1 )
all.equal( summary( b1n )[ -c( 3, 7, 20, 42 ) ], 
   summary( b1ni )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b1n, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b1ni, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b1n )
lrtest( b1n )
round( efficiencies( b1n ), 2 )
round( efficiencies( b1n, asInData = TRUE ), 2 )
all.equal( efficiencies( b1n ), efficiencies( b1ni ) )
all.equal( efficiencies( b1n, asInData = TRUE, minusU = FALSE ), 
   efficiencies( b1ni, asInData = TRUE, minusU = FALSE ) )
round( residuals( b1n ), 2 )
round( residuals( b1n, asInData = TRUE ), 2 )
all.equal( fitted( b1n, asInData = TRUE ) + residuals( b1n, asInData = TRUE ),
   log( naPanelData$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( naPanelData$PROD ) - fitted( b1n, asInData = TRUE ),
   residuals( b1n, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b1n )

## panel data with NA firms, error components frontier, truncNorm, timeEffect
b4n <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ),
   data = naPanelData, truncNorm = TRUE, timeEffect = TRUE )
b4ni <- sfa( log( PROD ) ~ ones + log( AREA ) + log( LABOR ) + log( NPK ) - 1,
   data = naPanelData, truncNorm = TRUE, timeEffect = TRUE )
all.equal( b4ni[ -c( 3, 7, 20, 42 ) ], b4n[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b4n, digits = 1 )
print( summary( b4n ), digits = 1 )
all.equal( summary( b4n )[ -c( 3, 7, 20, 42 ) ], 
   summary( b4ni )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b4n, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b4ni, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b4n )
lrtest( b4n )
round( efficiencies( b4n ), 2 )
round( efficiencies( b4n, asInData = TRUE ), 2 )
all.equal( efficiencies( b4n, minusU = FALSE ), 
   efficiencies( b4ni, minusU = FALSE ) )
all.equal( efficiencies( b4n, asInData = TRUE ), 
   efficiencies( b4ni, asInData = TRUE ) )
round( residuals( b4n ), 2 )
round( residuals( b4n, asInData = TRUE ), 2 )
all.equal( fitted( b4n, asInData = TRUE ) + residuals( b4n, asInData = TRUE ),
   log( naPanelData$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( naPanelData$PROD ) - fitted( b4n, asInData = TRUE ),
   residuals( b4n, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b4n )

## panel data with NA firms, efficiency effects frontier
b5n <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ) |
   EDYRS + BANRAT - 1, data = naPanelData )
b5ni <- sfa( log( PROD ) ~ ones + log( AREA ) + log( LABOR ) + log( NPK ) - 1 |
   EDYRS + BANRAT - 1, data = naPanelData )
all.equal( b5ni[ -c( 3, 7, 20, 42 ) ], b5n[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b5n, digits = 1 )
print( summary( b5n ), digits = 1 )
all.equal( summary( b5n )[ -c( 3, 7, 20, 42 ) ], 
   summary( b5ni )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b5n, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b5ni, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b5n )
lrtest( b5n )
printME( efficiencies( b5n, margEff = TRUE ) )
printME( efficiencies( b5n, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( b5n, margEff = TRUE ), 
   efficiencies( b5ni, margEff = TRUE ) )
all.equal( efficiencies( b5n, asInData = TRUE, margEff = TRUE, minusU = FALSE ), 
   efficiencies( b5ni, asInData = TRUE, margEff = TRUE, minusU = FALSE ) )
round( residuals( b5n ), 2 )
round( residuals( b5n, asInData = TRUE ), 2 )
all.equal( fitted( b5n, asInData = TRUE ) + residuals( b5n, asInData = TRUE ),
   log( naPanelData$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( naPanelData$PROD ) - fitted( b5n, asInData = TRUE ),
   residuals( b5n, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b5n )

## panel data with NA firms, efficiency effects frontier, zIntercept
b6n <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ) |
   EDYRS + BANRAT, data = naPanelData )
b6ni <- sfa( log( PROD ) ~ ones + log( AREA ) + log( LABOR ) + log( NPK ) - 1 |
      EDYRS + BANRAT, data = naPanelData )
all.equal( b6ni[ -c( 3, 7, 20, 42 ) ], b6n[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b6n, digits = 1 )
print( summary( b6n ), digits = 1 )
all.equal( summary( b6n )[ -c( 3, 7, 20, 42 ) ], 
   summary( b6ni )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b6n, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b6ni, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b6n )
lrtest( b6n )
printME( efficiencies( b6n, margEff = TRUE ) )
printME( efficiencies( b6n, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( b6n, margEff = TRUE, minusU = FALSE ), 
   efficiencies( b6ni, margEff = TRUE, minusU = FALSE ) )
all.equal( efficiencies( b6n, asInData = TRUE, margEff = TRUE ), 
   efficiencies( b6ni, asInData = TRUE, margEff = TRUE ) )
round( residuals( b6n ), 2 )
round( residuals( b6n, asInData = TRUE ), 2 )
all.equal( fitted( b6n, asInData = TRUE ) + residuals( b6n, asInData = TRUE ),
   log( naPanelData$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( naPanelData$PROD ) - fitted( b6n, asInData = TRUE ),
   residuals( b6n, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b6n )


## unbalanced panel data with time periods that have NAs for all firms
naTimePanelData <- riceProdPhilPanelUnb
naTimePanelData[ naTimePanelData$year == 2001, "PROD" ] <- NA
naTimePanelData[ naTimePanelData$year == 2004, "AREA" ] <- NA
naTimePanelData[ naTimePanelData$year == 1999, "EDYRS" ] <- NA

## panel data with NA years, error components frontier
b1t <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ),
   data = naTimePanelData )
b1ti <- sfa( log( PROD ) ~ ones + log( AREA ) + log( LABOR ) + log( NPK ) - 1,
   data = naTimePanelData )
all.equal( b1ti[ -c( 3, 7, 20, 42 ) ], b1t[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b1t, digits = 1 )
print( summary( b1t ), digits = 1 )
all.equal( summary( b1t )[ -c( 3, 7, 20, 42 ) ], 
   summary( b1ti )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b1t, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b1ti, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b1t )
lrtest( b1t )
round( efficiencies( b1t ), 2 )
round( efficiencies( b1t, asInData = TRUE ), 2 )
all.equal( efficiencies( b1t, minusU = FALSE ), 
   efficiencies( b1ti, minusU = FALSE ) )
all.equal( efficiencies( b1t, asInData = TRUE ), 
   efficiencies( b1ti, asInData = TRUE ) )
round( fitted( b1t ), 2 )
round( fitted( b1t, asInData = TRUE ), 2 )
round( residuals( b1t ), 2 )
round( residuals( b1t, asInData = TRUE ), 2 )
all.equal( fitted( b1t, asInData = TRUE ) + residuals( b1t, asInData = TRUE ),
   log( naTimePanelData$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( naTimePanelData$PROD ) - fitted( b1t, asInData = TRUE ),
   residuals( b1t, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b1t )

## panel data with NA years, error components frontier, truncNorm, timeEffect
b4t <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ),
   data = naTimePanelData, truncNorm = TRUE, timeEffect = TRUE )
b4ti <- sfa( log( PROD ) ~ ones + log( AREA ) + log( LABOR ) + log( NPK ) - 1,
   data = naTimePanelData, truncNorm = TRUE, timeEffect = TRUE )
all.equal( b4ti[ -c( 3, 7, 20, 42 ) ], b4t[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b4t, digits = 1 )
print( summary( b4t ), digits = 1 )
all.equal( summary( b4t )[ -c( 3, 7, 20, 42 ) ], 
   summary( b4ti )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b4t, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b4ti, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b4t )
lrtest( b4t )
round( efficiencies( b4t ), 2 )
round( efficiencies( b4t, asInData = TRUE ), 2 )
all.equal( efficiencies( b4t ), efficiencies( b4ti ) )
all.equal( efficiencies( b4t, asInData = TRUE, minusU = FALSE ), 
   efficiencies( b4ti, asInData = TRUE, minusU = FALSE ) )
round( residuals( b4t ), 2 )
round( residuals( b4t, asInData = TRUE ), 2 )
all.equal( fitted( b4t, asInData = TRUE ) + residuals( b4t, asInData = TRUE ),
   log( naTimePanelData$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( naTimePanelData$PROD ) - fitted( b4t, asInData = TRUE ),
   residuals( b4t, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b4t )

## panel data with NA years, efficiency effects frontier
b5t <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ) |
   EDYRS + BANRAT - 1, data = naTimePanelData )
b5ti <- sfa( log( PROD ) ~ ones + log( AREA ) + log( LABOR ) + log( NPK ) - 1 |
      EDYRS + BANRAT - 1, data = naTimePanelData )
all.equal( b5ti[ -c( 3, 7, 20, 42 ) ], b5t[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b5t, digits = 1 )
print( summary( b5t ), digits = 1 )
all.equal( summary( b5t )[ -c( 3, 7, 20, 42 ) ], 
   summary( b5ti )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b5t, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b5ti, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b5t )
lrtest( b5t )
printME( efficiencies( b5t, margEff = TRUE ) )
printME( efficiencies( b5t, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( b5t, margEff = TRUE, minusU = FALSE ), 
   efficiencies( b5ti, margEff = TRUE, minusU = FALSE ) )
all.equal( efficiencies( b5t, asInData = TRUE, margEff = TRUE ), 
   efficiencies( b5ti, asInData = TRUE, margEff = TRUE ) )
round( fitted( b5t ), 2 )
round( fitted( b5t, asInData = TRUE ), 2 )
round( residuals( b5t ), 2 )
round( residuals( b5t, asInData = TRUE ), 2 )
all.equal( fitted( b5t, asInData = TRUE ) + residuals( b5t, asInData = TRUE ),
   log( naTimePanelData$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( naTimePanelData$PROD ) - fitted( b5t, asInData = TRUE ),
   residuals( b5t, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b5t )

## panel data with NA years, efficiency effects frontier, zIntercept
b6t <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ) |
   EDYRS + BANRAT, data = naTimePanelData )
b6ti <- sfa( log( PROD ) ~ ones + log( AREA ) + log( LABOR ) + log( NPK ) - 1 |
      EDYRS + BANRAT, data = naTimePanelData )
all.equal( b6ti[ -c( 3, 7, 20, 42 ) ], b6t[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( b6t, digits = 1 )
print( summary( b6t ), digits = 1 )
all.equal( summary( b6t )[ -c( 3, 7, 20, 42 ) ], 
   summary( b6ti )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( b6t, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( b6ti, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
nobs( b6t )
lrtest( b6t )
printME( efficiencies( b6t, margEff = TRUE ) )
printME( efficiencies( b6t, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( b6t, margEff = TRUE ), 
   efficiencies( b6ti, margEff = TRUE ) )
all.equal( efficiencies( b6t, asInData = TRUE, margEff = TRUE, minusU = FALSE ), 
   efficiencies( b6ti, asInData = TRUE, margEff = TRUE, minusU = FALSE ) )
round( residuals( b6t ), 2 )
round( residuals( b6t, asInData = TRUE ), 2 )
all.equal( fitted( b6t, asInData = TRUE ) + residuals( b6t, asInData = TRUE ),
   log( naTimePanelData$PROD ), check.attributes = FALSE, tol = 1e-4 )
all.equal( log( naTimePanelData$PROD ) - fitted( b6t, asInData = TRUE ),
   residuals( b6t, asInData = TRUE ), check.attributes = FALSE, tol = 1e-4 )
printAll( b6t )


## translog frontiers
## cross-section data, error components frontier, translog
translog <- frontierQuad( data = front41Data, yName = "logOutput",
   xNames = c( "logCapital", "logLabour" ) )
print( translog, digits = 1 )
coef( translog, which = "start" )
round( coef( translog, which = "ols" ), 2 )
round( coef( translog, which = "grid" ), 2 )
round( coef( translog ), 2 )
round( coef( summary( translog ), which = "ols" ), 2 )
round( coef( summary( translog ) ), 2 )
round( vcov( translog ), 2 )
print( logLik( translog, which = "ols" ), digits = 4 )
print( logLik( translog ), digits = 4 )
nobs( translog )
print( summary( translog ), digits = 1 )
print( summary( translog, effMinusU = FALSE ), digits = 1 )
lrtest( translog )
round( efficiencies( translog ), 2 )
round( efficiencies( translog, asInData = TRUE ), 2 )
round( efficiencies( translog, minusU = FALSE ), 2 )
round( efficiencies( translog, asInData = TRUE, minusU = FALSE ), 2 )
round( fitted( translog ), 2 )
round( fitted( translog, asInData = TRUE ), 2 )
round( residuals( translog ), 2 )
round( residuals( translog, asInData = TRUE ), 2 )
all.equal( fitted( translog, asInData = TRUE ) + 
      residuals( translog, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )
translogEla <- elas( translog )
round( translogEla, 2 )
round( attributes( translogEla )$variance, 2 )
round( attributes( translogEla )$stdDev, 2 )
printAll( translog )

## cross-section data, error components frontier, translog, shifter
translogShift <- frontierQuad( yName = "logOutput",
   xNames = c( "logCapital", "logLabour" ), shifterNames = "firmNo",
   data = front41Data )
print( translogShift, digits = 1 )
coef( translogShift, which = "start" )
round( coef( translogShift, which = "ols" ), 2 )
round( coef( translogShift, which = "grid" ), 2 )
round( coef( translogShift ), 2 )
round( coef( summary( translogShift ), which = "ols" ), 2 )
round( coef( summary( translogShift ) ), 2 )
round( vcov( translogShift ), 2 )
print( logLik( translogShift, which = "ols" ), digits = 4 )
print( logLik( translogShift ), digits = 4 )
nobs( translogShift )
print( summary( translogShift ), digits = 1 )
lrtest( translogShift )
round( efficiencies( translogShift ), 2 )
round( efficiencies( translogShift, asInData = TRUE ), 2 )
round( residuals( translogShift ), 2 )
round( residuals( translogShift, asInData = TRUE ), 2 )
all.equal( fitted( translogShift, asInData = TRUE ) + 
      residuals( translogShift, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )
translogShiftEla <- elas( translogShift )
round( translogShiftEla, 2 )
round( attributes( translogShiftEla )$variance, 2 )
round( attributes( translogShiftEla )$stdDev, 2 )
printAll( translogShift )

## cross-section data, efficiency effects frontier, translog
translogZvar <- frontierQuad( yName = "logOutput",
   xNames = c( "logCapital", "logLabour" ), zNames = "firmNo",
   data = front41Data )
print( translogZvar, digits = 1 )
coef( translogZvar, which = "start" )
round( coef( translogZvar, which = "ols" ), 2 )
round( coef( translogZvar, which = "grid" ), 2 )
round( coef( translogZvar ), 2 )
round( coef( summary( translogZvar ), which = "ols" ), 2 )
round( coef( summary( translogZvar ) ), 2 )
round( vcov( translogZvar ), 2 )
print( logLik( translogZvar, which = "ols" ), digits = 4 )
print( logLik( translogZvar ), digits = 4 )
nobs( translogZvar )
print( summary( translogZvar ), digits = 1 )
print( summary( translogZvar, effMinusU = FALSE ), digits = 1 )
lrtest( translogZvar )
printME( efficiencies( translogZvar, margEff = TRUE ) )
printME( efficiencies( translogZvar, asInData = TRUE, margEff = TRUE ) )
printME( efficiencies( translogZvar, minusU = FALSE, margEff = TRUE ) )
printME( efficiencies( translogZvar, asInData = TRUE, minusU = FALSE, margEff = TRUE ) )
round( fitted( translogZvar ), 2 )
round( fitted( translogZvar, asInData = TRUE ), 2 )
round( residuals( translogZvar ), 2 )
round( residuals( translogZvar, asInData = TRUE ), 2 )
all.equal( fitted( translogZvar, asInData = TRUE ) + 
      residuals( translogZvar, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )
translogZvarEla <- elas( translogZvar )
round( translogZvarEla, 2 )
round( attributes( translogZvarEla )$variance, 2 )
round( attributes( translogZvarEla )$stdDev, 2 )
printAll( translogZvar )


#######  only an intercept as explanatory variable  #######
## cross-section data, error components frontier
oi1 <- sfa( logOutput ~ 1, data = front41Data )
oi1i <- sfa( logOutput ~ ones - 1, data = front41Data )
all.equal( oi1[ -c( 3, 7, 20, 42 ) ], oi1i[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( oi1 )
coef( oi1, which = "start" )
round( coef( oi1, which = "ols" ), 2 )
round( coef( oi1, which = "grid" ), 2 )
round( coef( oi1 ), 2 )
round( coef( summary( oi1 ), which = "ols" ), 2 )
round( coef( summary( oi1 ) ), 2 )
round( vcov( oi1 ), 2 )
print( logLik( oi1, which = "ols" ), digits = 4 )
print( logLik( oi1, which = "grid" ), digits = 4 )
print( logLik( oi1 ), digits = 4 )
nobs( oi1 )
print( summary( oi1 ), digits = 1 )
print( summary( oi1, effMinusU = FALSE ), digits = 1 )
all.equal( summary( oi1 )[ -c( 3, 7, 20, 42 ) ], 
   summary( oi1i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( oi1, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( oi1i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( oi1 )
round( efficiencies( oi1 ), 2 )
round( efficiencies( oi1, asInData = TRUE ), 2 )
all.equal( efficiencies( oi1 ), efficiencies( oi1i ) )
all.equal( efficiencies( oi1, asInData = TRUE, minusU = FALSE ), 
   efficiencies( oi1i, asInData = TRUE, minusU = FALSE ) )
round( fitted( oi1 ), 2 )
round( residuals( oi1 ), 2 )
all.equal( fitted( oi1, asInData = TRUE ) + residuals( oi1, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )
printAll( oi1 )

## cross-section data, error components frontier, truncNorm
oi2 <- sfa( logOutput ~ 1, data = front41Data, truncNorm = TRUE )
oi2i <- sfa( logOutput ~ ones - 1, data = front41Data, truncNorm = TRUE )
all.equal( oi2[ -c( 3, 7, 20, 42 ) ], oi2i[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( oi2 )
round( coef( oi2, which = "ols" ), 2 )
round( coef( oi2, which = "grid" ), 2 )
round( coef( oi2 ), 2 )
round( coef( summary( oi2 ), which = "ols" ), 2 )
round( coef( summary( oi2 ) ), 2 )
round( vcov( oi2 ), 2 )
print( logLik( oi2, which = "ols" ), digits = 4 )
print( logLik( oi2, which = "grid" ), digits = 4 )
print( logLik( oi2 ), digits = 4 )
nobs( oi2 )
print( summary( oi2 ), digits = 1 )
print( summary( oi2, effMinusU = FALSE ), digits = 1 )
all.equal( summary( oi2 )[ -c( 3, 7, 20, 42 ) ], 
   summary( oi2i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( oi2, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( oi2i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( oi2 )
round( efficiencies( oi2 ), 2 )
round( efficiencies( oi2, asInData = TRUE ), 2 )
all.equal( efficiencies( oi2 ), efficiencies( oi2i ) )
all.equal( efficiencies( oi2, asInData = TRUE, minusU = FALSE ), 
   efficiencies( oi2i, asInData = TRUE, minusU = FALSE ) )
round( fitted( oi2 ), 2 )
round( residuals( oi2 ), 2 )
all.equal( fitted( oi2, asInData = TRUE ) + residuals( oi2, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )
printAll( oi2 )

## cross-section data, efficiency effects frontier
oi5 <- sfa( logOutput ~ 1 | firmNo - 1, data = front41Data )
oi5i <- sfa( logOutput ~ ones - 1 | firmNo - 1, data = front41Data )
all.equal( oi5[ -c( 3, 7, 20, 42 ) ], oi5i[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( oi5 )
round( coef( oi5, which = "ols" ), 2 )
round( coef( oi5, which = "grid" ), 2 )
round( coef( oi5 ), 2 )
round( coef( summary( oi5 ), which = "ols" ), 2 )
round( coef( summary( oi5 ) ), 2 )
round( vcov( oi5 ), 2 )
print( logLik( oi5, which = "ols" ), digits = 4 )
print( logLik( oi5, which = "grid" ), digits = 4 )
print( logLik( oi5 ), digits = 4 )
nobs( oi5 )
print( summary( oi5 ), digits = 1 )
print( summary( oi5, effMinusU = FALSE ), digits = 1 )
all.equal( summary( oi5 )[ -c( 3, 7, 20, 42 ) ], 
   summary( oi5i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( oi5, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( oi5i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( oi5 )
printME( efficiencies( oi5, margEff = TRUE ) )
printME( efficiencies( oi5, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( oi5, minusU = FALSE, margEff = TRUE ), 
   efficiencies( oi5i, minusU = FALSE, margEff = TRUE ) )
all.equal( efficiencies( oi5, asInData = TRUE, margEff = TRUE ), 
   efficiencies( oi5i, asInData = TRUE, margEff = TRUE ) )
round( fitted( oi5 ), 2 )
round( residuals( oi5 ), 2 )
all.equal( fitted( oi5, asInData = TRUE ) + residuals( oi5, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )
printAll( oi5 )

## cross-section data, efficiency effects frontier, zIntercept
oi6 <- sfa( logOutput ~ 1 | firmNo, data = front41Data )
oi6i <- sfa( logOutput ~ ones - 1 | firmNo, data = front41Data )
all.equal( oi6[ -c( 3, 7, 20, 42 ) ], oi6i[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( oi6 )
round( coef( oi6, which = "ols" ), 2 )
round( coef( oi6, which = "grid" ), 2 )
round( coef( oi6 ), 2 )
round( coef( summary( oi6 ), which = "ols" ), 2 )
round( coef( summary( oi6 ) ), 2 )
round( vcov( oi6 ), 2 )
print( logLik( oi6, which = "ols" ), digits = 4 )
print( logLik( oi6, which = "grid" ), digits = 4 )
print( logLik( oi6 ), digits = 4 )
nobs( oi6 )
print( summary( oi6 ), digits = 1 )
print( summary( oi6, effMinusU = FALSE ), digits = 1 )
all.equal( summary( oi6 )[ -c( 3, 7, 20, 42 ) ], 
   summary( oi6i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( oi6, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( oi6i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( oi6 )
printME( efficiencies( oi6, margEff = TRUE ) )
printME( efficiencies( oi6, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( oi6, margEff = TRUE ), 
   efficiencies( oi6i, margEff = TRUE ) )
all.equal( efficiencies( oi6, asInData = TRUE, minusU = FALSE, margEff = TRUE ), 
   efficiencies( oi6i, asInData = TRUE, minusU = FALSE, margEff = TRUE ) )
round( fitted( oi6 ), 2 )
round( residuals( oi6 ), 2 )
all.equal( fitted( oi6, asInData = TRUE ) + residuals( oi6, asInData = TRUE ),
   front41Data$logOutput, check.attributes = FALSE, tol = 1e-4 )
printAll( oi6 )


######  only an intercept as explanatory variable: unbalanced panel data  ######
riceProdPhilPanelUnb$lProdNa <- riceProdPhilPanelUnb$lPROD
riceProdPhilPanelUnb$lProdNa[ riceProdPhilPanelUnb$farm == "F_28" ] <- NA
riceProdPhilPanelUnb$lProdNa[ riceProdPhilPanelUnb$year == 2002 ] <- NA
riceProdPhilPanelUnb$lProdNa[ (1:20) * 17 ] <- NA


## unbalanced panel data, error components frontier
oip1 <- sfa( lProdNa ~ 1, data = riceProdPhilPanelUnb )
oip1i <- sfa( lProdNa ~ ones - 1, data = riceProdPhilPanelUnb )
all.equal( oip1[ -c( 3, 7, 20, 42 ) ], oip1i[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( oip1 )
round( coef( oip1, which = "ols" ), 2 )
round( coef( oip1, which = "grid" ), 2 )
round( coef( oip1 ), 2 )
round( coef( summary( oip1 ), which = "ols" ), 2 )
round( coef( summary( oip1 ) ), 2 )
round( vcov( oip1 ), 2 )
print( logLik( oip1, which = "ols" ), digits = 4 )
print( logLik( oip1, which = "grid" ), digits = 4 )
print( logLik( oip1 ), digits = 4 )
nobs( oip1 )
print( summary( oip1 ), digits = 1 )
print( summary( oip1, effMinusU = FALSE ), digits = 1 )
all.equal( summary( oip1 )[ -c( 3, 7, 20, 42 ) ], 
   summary( oip1i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( oip1, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( oip1i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( oip1 )
round( efficiencies( oip1 ), 2 )
round( efficiencies( oip1, asInData = TRUE ), 2 )
all.equal( efficiencies( oip1 ), efficiencies( oip1i ) )
all.equal( efficiencies( oip1, asInData = TRUE, minusU = FALSE ), 
   efficiencies( oip1i, asInData = TRUE, minusU = FALSE ) )
round( fitted( oip1 ), 2 )
round( residuals( oip1 ), 2 )
all.equal( fitted( oip1, asInData = TRUE ) + residuals( oip1, asInData = TRUE ),
   riceProdPhilPanelUnb$lProdNa, check.attributes = FALSE, tol = 1e-4 )
printAll( oip1 )

## unbalanced panel data, error components frontier, truncNorm
oip2 <- sfa( lProdNa ~ 1, data = riceProdPhilPanelUnb, truncNorm = TRUE )
oip2i <- sfa( lProdNa ~ ones - 1, data = riceProdPhilPanelUnb, truncNorm = TRUE )
all.equal( oip2[ -c( 3, 7, 20, 42 ) ], oip2i[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( oip2 )
round( coef( oip2, which = "ols" ), 2 )
round( coef( oip2, which = "grid" ), 2 )
round( coef( oip2 ), 2 )
round( coef( summary( oip2 ), which = "ols" ), 2 )
round( coef( summary( oip2 ) ), 2 )
round( vcov( oip2 ), 2 )
print( logLik( oip2, which = "ols" ), digits = 4 )
print( logLik( oip2, which = "grid" ), digits = 4 )
print( logLik( oip2 ), digits = 4 )
nobs( oip2 )
print( summary( oip2 ), digits = 1 )
print( summary( oip2, effMinusU = FALSE ), digits = 1 )
all.equal( summary( oip2 )[ -c( 3, 7, 20, 42 ) ], 
   summary( oip2i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( oip2, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( oip2i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( oip2 )
round( efficiencies( oip2 ), 2 )
round( efficiencies( oip2, asInData = TRUE ), 2 )
all.equal( efficiencies( oip2 ), efficiencies( oip2i ) )
all.equal( efficiencies( oip2, asInData = TRUE, minusU = FALSE ), 
   efficiencies( oip2i, asInData = TRUE, minusU = FALSE ) )
round( fitted( oip2 ), 2 )
round( residuals( oip2 ), 2 )
all.equal( fitted( oip2, asInData = TRUE ) + residuals( oip2, asInData = TRUE ),
   riceProdPhilPanelUnb$lProdNa, check.attributes = FALSE, tol = 1e-4 )
printAll( oip2 )

## unbalanced panel data, error components frontier, timeEffect
oip3 <- sfa( lProdNa ~ 1, data = riceProdPhilPanelUnb, timeEffect = TRUE )
oip3i <- sfa( lProdNa ~ ones - 1, data = riceProdPhilPanelUnb, timeEffect = TRUE )
all.equal( oip3[ -c( 3, 7, 20, 42 ) ], oip3i[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( oip3 )
round( coef( oip3, which = "ols" ), 2 )
round( coef( oip3, which = "grid" ), 2 )
round( coef( oip3 ), 2 )
round( coef( summary( oip3 ), which = "ols" ), 2 )
round( coef( summary( oip3 ) ), 2 )
round( vcov( oip3 ), 2 )
print( logLik( oip3, which = "ols" ), digits = 4 )
print( logLik( oip3, which = "grid" ), digits = 4 )
print( logLik( oip3 ), digits = 4 )
nobs( oip3 )
print( summary( oip3 ), digits = 1 )
print( summary( oip3, effMinusU = FALSE ), digits = 1 )
all.equal( summary( oip3 )[ -c( 3, 7, 20, 42 ) ], 
   summary( oip3i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( oip3, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( oip3i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( oip3 )
round( efficiencies( oip3 ), 2 )
round( efficiencies( oip3, asInData = TRUE ), 2 )
all.equal( efficiencies( oip3 ), efficiencies( oip3i ) )
all.equal( efficiencies( oip3, asInData = TRUE, minusU = FALSE ), 
   efficiencies( oip3i, asInData = TRUE, minusU = FALSE ) )
round( fitted( oip3 ), 2 )
round( residuals( oip3 ), 2 )
all.equal( fitted( oip3, asInData = TRUE ) + residuals( oip3, asInData = TRUE ),
   riceProdPhilPanelUnb$lProdNa, check.attributes = FALSE, tol = 1e-4 )
printAll( oip3 )

## unbalanced panel data, error components frontier, truncNorm, timeEffect
oip4 <- sfa( lProdNa ~ 1, data = riceProdPhilPanelUnb, 
   truncNorm = TRUE, timeEffect = TRUE )
oip4i <- sfa( lProdNa ~ ones - 1, data = riceProdPhilPanelUnb, 
   truncNorm = TRUE, timeEffect = TRUE )
all.equal( oip4[ -c( 3, 7, 20, 42 ) ], oip4i[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( oip4 )
round( coef( oip4, which = "ols" ), 2 )
round( coef( oip4, which = "grid" ), 2 )
round( coef( oip4 ), 2 )
round( coef( summary( oip4 ), which = "ols" ), 2 )
round( coef( summary( oip4 ) ), 2 )
round( vcov( oip4 ), 2 )
print( logLik( oip4, which = "ols" ), digits = 4 )
print( logLik( oip4, which = "grid" ), digits = 4 )
print( logLik( oip4 ), digits = 4 )
nobs( oip4 )
print( summary( oip4 ), digits = 1 )
print( summary( oip4, effMinusU = FALSE ), digits = 1 )
all.equal( summary( oip4 )[ -c( 3, 7, 20, 42 ) ], 
   summary( oip4i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( oip4, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( oip4i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( oip4 )
round( efficiencies( oip4 ), 2 )
round( efficiencies( oip4, asInData = TRUE ), 2 )
all.equal( efficiencies( oip4 ), efficiencies( oip4i ) )
all.equal( efficiencies( oip4, asInData = TRUE, minusU = FALSE ), 
   efficiencies( oip4i, asInData = TRUE, minusU = FALSE ) )
round( fitted( oip4 ), 2 )
round( residuals( oip4 ), 2 )
all.equal( fitted( oip4, asInData = TRUE ) + residuals( oip4, asInData = TRUE ),
   riceProdPhilPanelUnb$lProdNa, check.attributes = FALSE, tol = 1e-4 )
printAll( oip4 )

## unbalanced panel data, efficiency effects frontier, zIntercept
oip5 <- sfa( lProdNa ~ 1 | EDYRS + BANRAT - 1, data = riceProdPhilPanelUnb )
oip5i <- sfa( lProdNa ~ ones - 1 | EDYRS + BANRAT - 1, data = riceProdPhilPanelUnb )
all.equal( oip5[ -c( 3, 7, 20, 42 ) ], oip5i[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( oip5 )
round( coef( oip5, which = "ols" ), 2 )
round( coef( oip5, which = "grid" ), 2 )
round( coef( oip5 ), 2 )
round( coef( summary( oip5 ), which = "ols" ), 2 )
round( coef( summary( oip5 ) ), 2 )
round( vcov( oip5 ), 2 )
print( logLik( oip5, which = "ols" ), digits = 4 )
print( logLik( oip5, which = "grid" ), digits = 4 )
print( logLik( oip5 ), digits = 4 )
nobs( oip5 )
print( summary( oip5 ), digits = 1 )
print( summary( oip5, effMinusU = FALSE ), digits = 1 )
all.equal( summary( oip5 )[ -c( 3, 7, 20, 42 ) ], 
   summary( oip5i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( oip5, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( oip5i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( oip5 )
printME( efficiencies( oip5, margEff = TRUE ) )
printME( efficiencies( oip5, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( oip5, margEff = TRUE ), 
   efficiencies( oip5i, margEff = TRUE ) )
all.equal( efficiencies( oip5, asInData = TRUE, minusU = FALSE, margEff = TRUE ), 
   efficiencies( oip5i, asInData = TRUE, minusU = FALSE, margEff = TRUE ) )
round( fitted( oip5 ), 2 )
round( residuals( oip5 ), 2 )
all.equal( fitted( oip5, asInData = TRUE ) + residuals( oip5, asInData = TRUE ),
   riceProdPhilPanelUnb$lProdNa, check.attributes = FALSE, tol = 1e-4 )
printAll( oip5 )

## unbalanced panel data, efficiency effects frontier, zIntercept
oip6 <- sfa( lProdNa ~ 1 | EDYRS + BANRAT, data = riceProdPhilPanelUnb )
oip6i <- sfa( lProdNa ~ ones - 1 | EDYRS + BANRAT, data = riceProdPhilPanelUnb )
all.equal( oip6[ -c( 3, 7, 20, 42 ) ], oip6i[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE, tol = 1e-4 )
print( oip6 )
round( coef( oip6, which = "ols" ), 2 )
round( coef( oip6, which = "grid" ), 2 )
round( coef( oip6 ), 2 )
round( coef( summary( oip6 ), which = "ols" ), 2 )
round( coef( summary( oip6 ) ), 2 )
round( vcov( oip6 ), 2 )
print( logLik( oip6, which = "ols" ), digits = 4 )
print( logLik( oip6, which = "grid" ), digits = 4 )
print( logLik( oip6 ), digits = 4 )
nobs( oip6 )
print( summary( oip6 ), digits = 1 )
print( summary( oip6, effMinusU = FALSE ), digits = 1 )
all.equal( summary( oip6 )[ -c( 3, 7, 20, 42 ) ], 
   summary( oip6i )[ -c( 3, 7, 20, 42 ) ], check.attributes = FALSE )
all.equal( summary( oip6, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   summary( oip6i, effMinusU = FALSE )[ -c( 3, 7, 20, 42 ) ], 
   check.attributes = FALSE )
lrtest( oip6 )
printME( efficiencies( oip6, margEff = TRUE ) )
printME( efficiencies( oip6, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( oip6, margEff = TRUE ), 
   efficiencies( oip6i, margEff = TRUE ) )
all.equal( efficiencies( oip6, asInData = TRUE, minusU = FALSE, margEff = TRUE ), 
   efficiencies( oip6i, asInData = TRUE, minusU = FALSE, margEff = TRUE ) )
round( fitted( oip6 ), 2 )
round( residuals( oip6 ), 2 )
all.equal( fitted( oip6, asInData = TRUE ) + residuals( oip6, asInData = TRUE ),
   riceProdPhilPanelUnb$lProdNa, check.attributes = FALSE, tol = 1e-4 )
printAll( oip6 )


#######  no explanatory variable  #######
front41Data$ynx <- front41Data$logOutput - mean( front41Data$logOutput ) - 0.1

## cross-section data, error components frontier
nx1 <- sfa( ynx ~ -1, data = front41Data )
print( nx1 )
coef( nx1, which = "start" )
round( coef( nx1, which = "ols" ), 2 )
round( coef( nx1, which = "grid" ), 2 )
round( coef( nx1 ), 2 )
round( coef( summary( nx1 ), which = "ols" ), 2 )
round( coef( summary( nx1 ) ), 2 )
round( vcov( nx1 ), 2 )
print( logLik( nx1, which = "ols" ), digits = 4 )
print( logLik( nx1, which = "grid" ), digits = 4 )
print( logLik( nx1 ), digits = 4 )
nobs( nx1 )
print( summary( nx1 ), digits = 1 )
print( summary( nx1, effMinusU = FALSE ), digits = 1 )
lrtest( nx1 )
round( efficiencies( nx1 ), 2 )
round( efficiencies( nx1, asInData = TRUE ), 2 )
round( fitted( nx1 ), 2 )
round( residuals( nx1 ), 2 )
all.equal( fitted( nx1, asInData = TRUE ) + residuals( nx1, asInData = TRUE ),
   front41Data$ynx, check.attributes = FALSE, tol = 1e-4 )
printAll( nx1 )

## cross-section data, error components frontier, truncNorm
nx2 <- sfa( ynx ~ -1, data = front41Data, truncNorm = TRUE )
print( nx2 )
round( coef( nx2, which = "ols" ), 2 )
round( coef( nx2, which = "grid" ), 2 )
round( coef( nx2 ), 2 )
round( coef( summary( nx2 ), which = "ols" ), 2 )
round( coef( summary( nx2 ) ), 2 )
round( vcov( nx2 ), 2 )
print( logLik( nx2, which = "ols" ), digits = 4 )
print( logLik( nx2, which = "grid" ), digits = 4 )
print( logLik( nx2 ), digits = 4 )
nobs( nx2 )
print( summary( nx2 ), digits = 1 )
print( summary( nx2, effMinusU = FALSE ), digits = 1 )
lrtest( nx2 )
round( efficiencies( nx2 ), 2 )
round( efficiencies( nx2, asInData = TRUE ), 2 )
round( fitted( nx2 ), 2 )
round( residuals( nx2 ), 2 )
all.equal( fitted( nx2, asInData = TRUE ) + residuals( nx2, asInData = TRUE ),
   front41Data$ynx, check.attributes = FALSE, tol = 1e-4 )
printAll( nx2 )

## cross-section data, efficiency effects frontier
nx5 <- sfa( ynx ~ -1 | firmNo - 1, data = front41Data )
print( nx5 )
round( coef( nx5, which = "ols" ), 2 )
round( coef( nx5, which = "grid" ), 2 )
round( coef( nx5 ), 2 )
round( coef( summary( nx5 ), which = "ols" ), 2 )
round( coef( summary( nx5 ) ), 2 )
round( vcov( nx5 ), 2 )
print( logLik( nx5, which = "ols" ), digits = 4 )
print( logLik( nx5, which = "grid" ), digits = 4 )
print( logLik( nx5 ), digits = 4 )
nobs( nx5 )
print( summary( nx5 ), digits = 1 )
print( summary( nx5, effMinusU = FALSE ), digits = 1 )
lrtest( nx5 )
printME( efficiencies( nx5, margEff = TRUE ) )
printME( efficiencies( nx5, asInData = TRUE, margEff = TRUE ) )
round( fitted( nx5 ), 2 )
round( residuals( nx5 ), 2 )
all.equal( fitted( nx5, asInData = TRUE ) + residuals( nx5, asInData = TRUE ),
   front41Data$ynx, check.attributes = FALSE, tol = 1e-4 )
printAll( nx5 )

## cross-section data, efficiency effects frontier, zIntercept
nx6 <- sfa( ynx ~ -1 | firmNo, data = front41Data )
print( nx6 )
round( coef( nx6, which = "ols" ), 2 )
round( coef( nx6, which = "grid" ), 2 )
round( coef( nx6 ), 2 )
round( coef( summary( nx6 ), which = "ols" ), 2 )
round( coef( summary( nx6 ) ), 2 )
round( vcov( nx6 ), 2 )
print( logLik( nx6, which = "ols" ), digits = 4 )
print( logLik( nx6, which = "grid" ), digits = 4 )
print( logLik( nx6 ), digits = 4 )
nobs( nx6 )
print( summary( nx6 ), digits = 1 )
print( summary( nx6, effMinusU = FALSE ), digits = 1 )
lrtest( nx6 )
printME( efficiencies( nx6, margEff = TRUE ) )
printME( efficiencies( nx6, asInData = TRUE, margEff = TRUE ) )
round( fitted( nx6 ), 2 )
round( residuals( nx6 ), 2 )
all.equal( fitted( nx6, asInData = TRUE ) + residuals( nx6, asInData = TRUE ),
   front41Data$ynx, check.attributes = FALSE, tol = 1e-4 )
printAll( nx6 )


#######  no explanatory variable  #######
riceProdPhilPanelUnb$ynx <- riceProdPhilPanelUnb$lProdNa - 
   mean( riceProdPhilPanelUnb$lProdNa, na.rm = TRUE ) - 0.1

## unbalanced panel data, error components frontier
nxp1 <- sfa( ynx ~ -1, data = riceProdPhilPanelUnb )
nxp1r <- sfa( ynx ~ -1,
   data = riceProdPhilPanelUnb[ !is.na( riceProdPhilPanelUnb$ynx ), ] )
all.equal( nxp1[ -c(42, 43) ], nxp1r[ -c(42, 43) ] )
print( nxp1 )
round( coef( nxp1, which = "ols" ), 2 )
round( coef( nxp1, which = "grid" ), 2 )
round( coef( nxp1 ), 2 )
round( coef( summary( nxp1 ), which = "ols" ), 2 )
round( coef( summary( nxp1 ) ), 2 )
round( vcov( nxp1 ), 2 )
print( logLik( nxp1, which = "ols" ), digits = 4 )
print( logLik( nxp1, which = "grid" ), digits = 4 )
print( logLik( nxp1 ), digits = 4 )
nobs( nxp1 )
print( summary( nxp1 ), digits = 1 )
print( summary( nxp1, effMinusU = FALSE ), digits = 1 )
lrtest( nxp1 )
round( efficiencies( nxp1 ), 2 )
round( efficiencies( nxp1, asInData = TRUE ), 2 )
all.equal( efficiencies( nxp1r, asInData = TRUE ),
   efficiencies( nxp1, asInData = TRUE )[ !is.na( riceProdPhilPanelUnb$ynx ) ] )
round( fitted( nxp1 ), 2 )
round( residuals( nxp1 ), 2 )
all.equal( fitted( nxp1, asInData = TRUE ) + residuals( nxp1, asInData = TRUE ),
   riceProdPhilPanelUnb$ynx, check.attributes = FALSE, tol = 1e-4 )
all.equal( fitted( nxp1r, asInData = TRUE ) + residuals( nxp1r, asInData = TRUE ),
   riceProdPhilPanelUnb$ynx[ !is.na( riceProdPhilPanelUnb$ynx ) ], 
   check.attributes = FALSE, tol = 1e-4 )
printAll( nxp1 )

## unbalanced panel data, error components frontier, truncNorm
nxp2 <- sfa( ynx ~ -1, data = riceProdPhilPanelUnb, truncNorm = TRUE )
nxp2r <- sfa( ynx ~ -1,
   data = riceProdPhilPanelUnb[ !is.na( riceProdPhilPanelUnb$ynx ), ], 
   truncNorm = TRUE )
all.equal( nxp2[ -c(42, 43) ], nxp2r[ -c(42, 43) ] )
print( nxp2 )
round( coef( nxp2, which = "ols" ), 2 )
round( coef( nxp2, which = "grid" ), 2 )
round( coef( nxp2 ), 2 )
round( coef( summary( nxp2 ), which = "ols" ), 2 )
round( coef( summary( nxp2 ) ), 2 )
round( vcov( nxp2 ), 2 )
print( logLik( nxp2, which = "ols" ), digits = 4 )
print( logLik( nxp2, which = "grid" ), digits = 4 )
print( logLik( nxp2 ), digits = 4 )
nobs( nxp2 )
print( summary( nxp2 ), digits = 1 )
print( summary( nxp2, effMinusU = FALSE ), digits = 1 )
lrtest( nxp2 )
round( efficiencies( nxp2 ), 2 )
round( efficiencies( nxp2, asInData = TRUE ), 2 )
all.equal( efficiencies( nxp2r, asInData = TRUE ),
   efficiencies( nxp2, asInData = TRUE )[ !is.na( riceProdPhilPanelUnb$ynx ) ] )
round( fitted( nxp2 ), 2 )
round( residuals( nxp2 ), 2 )
all.equal( fitted( nxp2, asInData = TRUE ) + residuals( nxp2, asInData = TRUE ),
   riceProdPhilPanelUnb$ynx, check.attributes = FALSE, tol = 1e-4 )
all.equal( fitted( nxp2r, asInData = TRUE ) + residuals( nxp2r, asInData = TRUE ),
   riceProdPhilPanelUnb$ynx[ !is.na( riceProdPhilPanelUnb$ynx ) ], 
   check.attributes = FALSE, tol = 1e-4 )
printAll( nxp2 )

## unbalanced panel data, error components frontier, timeEffect
nxp3 <- sfa( ynx ~ -1, data = riceProdPhilPanelUnb, timeEffect = TRUE )
nxp3r <- sfa( ynx ~ -1,
   data = riceProdPhilPanelUnb[ !is.na( riceProdPhilPanelUnb$ynx ), ], 
   timeEffect = TRUE )
all.equal( nxp3[ -c(42, 43) ], nxp3r[ -c(42, 43) ] )
print( nxp3 )
round( coef( nxp3, which = "ols" ), 2 )
round( coef( nxp3, which = "grid" ), 2 )
round( coef( nxp3 ), 2 )
round( coef( summary( nxp3 ), which = "ols" ), 2 )
round( coef( summary( nxp3 ) ), 2 )
round( vcov( nxp3 ), 2 )
print( logLik( nxp3, which = "ols" ), digits = 4 )
print( logLik( nxp3, which = "grid" ), digits = 4 )
print( logLik( nxp3 ), digits = 4 )
nobs( nxp3 )
print( summary( nxp3 ), digits = 1 )
print( summary( nxp3, effMinusU = FALSE ), digits = 1 )
lrtest( nxp3 )
round( efficiencies( nxp3 ), 2 )
round( efficiencies( nxp3, asInData = TRUE ), 2 )
all.equal( efficiencies( nxp3r, asInData = TRUE ),
   efficiencies( nxp3, asInData = TRUE )[ !is.na( riceProdPhilPanelUnb$ynx ) ] )
round( fitted( nxp3 ), 2 )
round( residuals( nxp3 ), 2 )
all.equal( fitted( nxp3, asInData = TRUE ) + residuals( nxp3, asInData = TRUE ),
   riceProdPhilPanelUnb$ynx, check.attributes = FALSE, tol = 1e-4 )
all.equal( fitted( nxp3r, asInData = TRUE ) + residuals( nxp3r, asInData = TRUE ),
   riceProdPhilPanelUnb$ynx[ !is.na( riceProdPhilPanelUnb$ynx ) ], 
   check.attributes = FALSE, tol = 1e-4 )
printAll( nxp3 )

## unbalanced panel data, error components frontier, truncNorm, timeEffect
nxp4 <- sfa( ynx ~ -1, data = riceProdPhilPanelUnb, 
   truncNorm = TRUE, timeEffect = TRUE )
nxp4r <- sfa( ynx ~ -1,
   data = riceProdPhilPanelUnb[ !is.na( riceProdPhilPanelUnb$ynx ), ], 
   truncNorm = TRUE, timeEffect = TRUE )
all.equal( nxp4[ -c(42, 43) ], nxp4r[ -c(42, 43) ] )
print( nxp4 )
round( coef( nxp4, which = "ols" ), 2 )
round( coef( nxp4, which = "grid" ), 2 )
round( coef( nxp4 ), 2 )
round( coef( summary( nxp4 ), which = "ols" ), 2 )
round( coef( summary( nxp4 ) ), 2 )
round( vcov( nxp4 ), 2 )
print( logLik( nxp4, which = "ols" ), digits = 4 )
print( logLik( nxp4, which = "grid" ), digits = 4 )
print( logLik( nxp4 ), digits = 4 )
nobs( nxp4 )
print( summary( nxp4 ), digits = 1 )
print( summary( nxp4, effMinusU = FALSE ), digits = 1 )
lrtest( nxp4 )
round( efficiencies( nxp4 ), 2 )
round( efficiencies( nxp4, asInData = TRUE ), 2 )
all.equal( efficiencies( nxp4r, asInData = TRUE ),
   efficiencies( nxp4, asInData = TRUE )[ !is.na( riceProdPhilPanelUnb$ynx ) ] )
round( fitted( nxp4 ), 2 )
round( residuals( nxp4 ), 2 )
all.equal( fitted( nxp4, asInData = TRUE ) + residuals( nxp4, asInData = TRUE ),
   riceProdPhilPanelUnb$ynx, check.attributes = FALSE, tol = 1e-4 )
all.equal( fitted( nxp4r, asInData = TRUE ) + residuals( nxp4r, asInData = TRUE ),
   riceProdPhilPanelUnb$ynx[ !is.na( riceProdPhilPanelUnb$ynx ) ], 
   check.attributes = FALSE, tol = 1e-4 )
printAll( nxp4 )

## unbalanced panel data, efficiency effects frontier, zIntercept
nxp5 <- sfa( ynx ~ -1 | EDYRS + BANRAT - 1, data = riceProdPhilPanelUnb )
nxp5r <- sfa( ynx ~ -1 | EDYRS + BANRAT - 1,
   data = riceProdPhilPanelUnb[ !is.na( riceProdPhilPanelUnb$ynx ), ] )
all.equal( nxp5[ -c(42, 43) ], nxp5r[ -c(42, 43) ] )
print( nxp5 )
round( coef( nxp5, which = "ols" ), 2 )
round( coef( nxp5, which = "grid" ), 2 )
round( coef( nxp5 ), 2 )
round( coef( summary( nxp5 ), which = "ols" ), 2 )
round( coef( summary( nxp5 ) ), 2 )
round( vcov( nxp5 ), 2 )
print( logLik( nxp5, which = "ols" ), digits = 4 )
print( logLik( nxp5, which = "grid" ), digits = 4 )
print( logLik( nxp5 ), digits = 4 )
nobs( nxp5 )
print( summary( nxp5 ), digits = 1 )
print( summary( nxp5, effMinusU = FALSE ), digits = 1 )
lrtest( nxp5 )
printME( efficiencies( nxp5, margEff = TRUE ) )
printME( efficiencies( nxp5, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( nxp5r, asInData = TRUE ),
   efficiencies( nxp5, asInData = TRUE )[ !is.na( riceProdPhilPanelUnb$ynx ) ] )
all.equal( 
   attr( efficiencies( nxp5r, asInData = TRUE, margEff = TRUE ), "margEff" ),
   attr( efficiencies( nxp5, asInData = TRUE, margEff = TRUE ), "margEff" )[ 
      !is.na( riceProdPhilPanelUnb$ynx ), ] )
round( fitted( nxp5 ), 2 )
round( residuals( nxp5 ), 2 )
all.equal( fitted( nxp5, asInData = TRUE ) + residuals( nxp5, asInData = TRUE ),
   riceProdPhilPanelUnb$ynx, check.attributes = FALSE, tol = 1e-4 )
all.equal( fitted( nxp5r, asInData = TRUE ) + residuals( nxp5r, asInData = TRUE ),
   riceProdPhilPanelUnb$ynx[ !is.na( riceProdPhilPanelUnb$ynx ) ], 
   check.attributes = FALSE, tol = 1e-4 )
printAll( nxp5 )

## unbalanced panel data, efficiency effects frontier, zIntercept
nxp6 <- sfa( ynx ~ -1 | EDYRS + BANRAT, data = riceProdPhilPanelUnb )
nxp6r <- sfa( ynx ~ -1 | EDYRS + BANRAT,
   data = riceProdPhilPanelUnb[ !is.na( riceProdPhilPanelUnb$ynx ), ] )
all.equal( nxp6[ -c(42, 43) ], nxp6r[ -c(42, 43) ] )
print( nxp6 )
round( coef( nxp6, which = "ols" ), 2 )
round( coef( nxp6, which = "grid" ), 2 )
round( coef( nxp6 ), 2 )
round( coef( summary( nxp6 ), which = "ols" ), 2 )
round( coef( summary( nxp6 ) ), 2 )
round( vcov( nxp6 ), 2 )
print( logLik( nxp6, which = "ols" ), digits = 4 )
print( logLik( nxp6, which = "grid" ), digits = 4 )
print( logLik( nxp6 ), digits = 4 )
nobs( nxp6 )
print( summary( nxp6 ), digits = 1 )
print( summary( nxp6, effMinusU = FALSE ), digits = 1 )
lrtest( nxp6 )
printME( efficiencies( nxp6, margEff = TRUE ) )
printME( efficiencies( nxp6, asInData = TRUE, margEff = TRUE ) )
all.equal( efficiencies( nxp6r, asInData = TRUE ),
   efficiencies( nxp6, asInData = TRUE )[ !is.na( riceProdPhilPanelUnb$ynx ) ] )
all.equal( 
   attr( efficiencies( nxp6r, asInData = TRUE, margEff = TRUE ), "margEff" ),
   attr( efficiencies( nxp6, asInData = TRUE, margEff = TRUE ), "margEff" )[ 
      !is.na( riceProdPhilPanelUnb$ynx ), ] )
round( fitted( nxp6 ), 2 )
round( residuals( nxp6 ), 2 )
all.equal( fitted( nxp6, asInData = TRUE ) + residuals( nxp6, asInData = TRUE ),
   riceProdPhilPanelUnb$ynx, check.attributes = FALSE, tol = 1e-4 )
all.equal( fitted( nxp6r, asInData = TRUE ) + residuals( nxp6r, asInData = TRUE ),
   riceProdPhilPanelUnb$ynx[ !is.na( riceProdPhilPanelUnb$ynx ) ], 
   check.attributes = FALSE, tol = 1e-4 )
printAll( nxp6 )


###  manual intercept and changing the values of the constant explanatory variable  ###
front41Data$tens <- 10
sa1iTens <- sfa( logOutput ~ tens + logCapital + logLabour - 1, 
   data = front41Data )
round( sa1iTens$gridAdj, 3 )
all.equal( coef( sa1iTens, which = "ols" ) * c( 10, 1, 1, 1 ), 
   coef( sa1i, which = "ols" ), check.attributes = FALSE )
all.equal( coef( sa1iTens, which = "grid" ) * c( 10, 1, 1, 1, 1 ), 
   coef( sa1i, which = "grid" ), check.attributes = FALSE )
all.equal( coef( sa1iTens ) * c( 10, 1, 1, 1, 1 ), coef( sa1i ), 
   check.attributes = FALSE, tol = 1e-4 )

front41Data$hundreds <- 100
sa1iHundreds <- sfa( logOutput ~ hundreds + logCapital + logLabour - 1, 
   data = front41Data )
round( sa1iHundreds$gridAdj, 3 )
all.equal( coef( sa1iHundreds, which = "ols" ) * c( 100, 1, 1, 1 ), 
   coef( sa1i, which = "ols" ), check.attributes = FALSE )
all.equal( coef( sa1iHundreds, which = "grid" ) * c( 100, 1, 1, 1, 1 ), 
   coef( sa1i, which = "grid" ), check.attributes = FALSE )
all.equal( coef( sa1iHundreds ) * c( 100, 1, 1, 1, 1 ), coef( sa1i ), 
   check.attributes = FALSE, tol = 1e-4 )


###  manual intercept and changing the order of the explanatory variables  ###
sa1iCap <- sfa( logOutput ~ logCapital + logLabour + ones - 1, data = front41Data )
round( sa1iCap$gridAdj, 3 )
all.equal( coef( sa1iCap, which = "ols" )[ c(3,1,2,4) ], 
   coef( sa1i, which = "ols" ) )
all.equal( logLik( sa1iCap, which = "ols" ), logLik( sa1i, which = "ols" ) )
all.equal( coef( sa1iCap, which = "grid" )[ c(3,1,2,4,5) ],
   coef( sa1i, which = "grid" ), check.attributes = FALSE )
all.equal( logLik( sa1iCap, which = "grid" ), logLik( sa1i, which = "grid" ) )
all.equal( coef( sa1iCap )[ c(3,1,2,4,5) ], coef( sa1i ), tol = 1e-5 )
all.equal( logLik( sa1iCap ), logLik( sa1i ) )

sa1iLab <- sfa( logOutput ~ logLabour + logCapital + ones - 1, data = front41Data )
round( sa1iLab$gridAdj, 3 )
all.equal( coef( sa1iLab, which = "ols" )[ c(3,2,1,4) ], 
   coef( sa1i, which = "ols" ) )
all.equal( logLik( sa1iLab, which = "ols" ), logLik( sa1i, which = "ols" ) )
all.equal( coef( sa1iLab, which = "grid" )[ c(3,2,1,4,5) ],
   coef( sa1i, which = "grid" ) )
all.equal( logLik( sa1iCap, which = "grid" ), logLik( sa1i, which = "grid" ) )
all.equal( coef( sa1iLab )[ c(3,2,1,4,5) ], coef( sa1i ), tol = 1e-5 )
all.equal( logLik( sa1iLab ), logLik( sa1i ) )

sa5iCap <- sfa( logOutput ~ logCapital + logLabour + ones - 1 | firmNo - 1, 
   data = front41Data )
round( sa5iCap$gridAdj, 3 )
all.equal( coef( sa5iCap, which = "ols" ), coef( sa1iCap, which = "ols" ) )
all.equal( coef( sa5iCap, which = "grid" ), coef( sa1iCap, which = "grid" ) )
all.equal( coef( sa5iCap )[ c(3,1,2,4,5,6) ], coef( saa1i ), tol = 1e-5 )

sa5iLab <- sfa( logOutput ~ logLabour + logCapital + ones - 1 | firmNo - 1, 
   data = front41Data )
round( sa5iLab$gridAdj, 3 )
all.equal( coef( sa5iLab, which = "ols" ), coef( sa1iLab, which = "ols" ) )
all.equal( coef( sa5iLab, which = "grid" ), coef( sa1iLab, which = "grid" ) )
all.equal( coef( sa5iLab )[ c(3,2,1,4,5,6) ], coef( saa1i ), tol = 1e-5 )

###  no intercept at all and changing the order of the explanatory variables  ###
front41Data$mlogOutput <- front41Data$logOutput - mean( front41Data$logOutput )

sa1iLabCap <- sfa( mlogOutput ~ logLabour + logCapital - 1, data = front41Data )
summary( sa1iLabCap )
round( sa1iLabCap$gridAdj, 3 )
coef( sa1iLabCap, which = "ols" )
logLik( sa1iLabCap, which = "ols" )
coef( sa1iLabCap, which = "grid" )
logLik( sa1iLabCap, which = "grid" )
coef( sa1iLabCap )
logLik( sa1iLabCap  )

sa1iCapLab <- sfa( mlogOutput ~ logCapital + logLabour - 1, data = front41Data )
summary( sa1iCapLab )
round( sa1iCapLab$gridAdj, 3 )
coef( sa1iCapLab, which = "ols" )
all.equal( coef( sa1iCapLab, which = "ols" )[ c(2,1,3) ], 
   coef( sa1iLabCap, which = "ols" ) )
logLik( sa1iCapLab, which = "ols" )
all.equal( logLik( sa1iCapLab, which = "ols" ), 
   logLik( sa1iCapLab, which = "ols" ) )
coef( sa1iCapLab, which = "grid" )
all.equal( coef( sa1iCapLab, which = "grid" )[ c(2,1,3,4) ],
   coef( sa1iLabCap, which = "grid" ) )
logLik( sa1iCapLab, which = "grid" )
all.equal( logLik( sa1iCapLab, which = "grid" ), 
   logLik( sa1iLabCap, which = "grid" ) )
coef( sa1iCapLab )
all.equal( coef( sa1iCapLab )[ c(2,1,3,4) ], coef( sa1iLabCap ) )
logLik( sa1iCapLab )
all.equal( logLik( sa1iCapLab ), logLik( sa1iLabCap ) )

sa5iCapLab <- sfa( mlogOutput ~ logCapital + logLabour - 1 | firmNo - 1, 
   data = front41Data )
round( sa5iCapLab$gridAdj, 3 )
all.equal( coef( sa5iCapLab, which = "ols" ), coef( sa1iCapLab, which = "ols" ) )
all.equal( coef( sa5iCapLab, which = "grid" ), coef( sa1iCapLab, which = "grid" ) )
coef( sa5iCapLab )

sa5iLabCap <- sfa( mlogOutput ~ logLabour + logCapital - 1 | firmNo - 1, 
   data = front41Data )
round( sa5iLabCap$gridAdj, 3 )
all.equal( coef( sa5iLabCap, which = "ols" ), coef( sa1iLabCap, which = "ols" ) )
all.equal( coef( sa5iLabCap, which = "grid" ), coef( sa1iLabCap, which = "grid" ) )
all.equal( coef( sa5iLabCap )[ c(2,1,3,4,5) ], coef( sa5iCapLab ), tol = 1e-5 )


###  no intercept at all and changing the scale of the explanatory variables  ###
front41Data$logCapitalTen <- 10 * front41Data$logCapital
front41Data$logLabourTen <- 10 * front41Data$logLabour

sa1iLabTenCap <- sfa( mlogOutput ~ logLabourTen + logCapital - 1, 
   data = front41Data )
round( sa1iLabTenCap$gridAdj, 3 )
all.equal( coef( sa1iLabTenCap, which = "ols" ) * c( 10, 1, 1 ), 
   coef( sa1iLabCap, which = "ols" ), check.attributes = FALSE )
all.equal( coef( sa1iLabTenCap, which = "grid" ) * c( 10, 1, 1, 1 ), 
   coef( sa1iLabCap, which = "grid" ), check.attributes = FALSE )
all.equal( coef( sa1iLabTenCap ) * c( 10, 1, 1, 1 ), coef( sa1iLabCap ),
   check.attributes = FALSE, tol = 1e-5)

sa1iLabCapTen <- sfa( mlogOutput ~ logLabour + logCapitalTen - 1, 
   data = front41Data )
round( sa1iLabCapTen$gridAdj, 3 )
all.equal( coef( sa1iLabCapTen, which = "ols" ) * c( 1, 10, 1 ), 
   coef( sa1iLabCap, which = "ols" ), check.attributes = FALSE )
all.equal( coef( sa1iLabCapTen, which = "grid" ) * c( 1, 10, 1, 1 ), 
   coef( sa1iLabCap, which = "grid" ), check.attributes = FALSE )
all.equal( coef( sa1iLabCapTen ) * c( 1, 10, 1, 1 ), coef( sa1iLabCap ),
   check.attributes = FALSE, tol = 1e-5 )

sa1iLabTenCapTen <- sfa( mlogOutput ~ logLabourTen + logCapitalTen - 1, 
   data = front41Data )
round( sa1iLabTenCapTen$gridAdj, 3 )
all.equal( coef( sa1iLabTenCapTen, which = "ols" ) * c( 10, 10, 1 ), 
   coef( sa1iLabCap, which = "ols" ), check.attributes = FALSE )
all.equal( coef( sa1iLabTenCapTen, which = "grid" ) * c( 1, 10, 1, 1 ), 
   coef( sa1iLabTenCap, which = "grid" ), check.attributes = FALSE )
all.equal( coef( sa1iLabTenCapTen ) * c( 1, 10, 1, 1 ), coef( sa1iLabTenCap ),
   check.attributes = FALSE, tol = 1e-4 )

sa1iCapTenLab <- sfa( mlogOutput ~ logCapitalTen + logLabour - 1, 
   data = front41Data )
round( sa1iCapTenLab$gridAdj, 3 )
all.equal( coef( sa1iCapTenLab, which = "ols" ) * c( 10, 1, 1 ), 
   coef( sa1iCapLab, which = "ols" ), check.attributes = FALSE )
all.equal( coef( sa1iCapTenLab, which = "grid" ) * c( 10, 1, 1, 1 ), 
   coef( sa1iCapLab, which = "grid" ), check.attributes = FALSE )
all.equal( coef( sa1iCapTenLab ) * c( 10, 1, 1, 1 ), coef( sa1iCapLab ),
   check.attributes = FALSE, tol = 1e-4 )

sa1iCapLabTen <- sfa( mlogOutput ~ logCapital + logLabourTen - 1, 
   data = front41Data )
round( sa1iCapLabTen$gridAdj, 3 )
all.equal( coef( sa1iCapLabTen, which = "ols" ) * c( 1, 10, 1 ), 
   coef( sa1iCapLab, which = "ols" ), check.attributes = FALSE )
all.equal( coef( sa1iCapLabTen, which = "grid" ) * c( 1, 10, 1, 1 ), 
   coef( sa1iCapLab, which = "grid" ), check.attributes = FALSE )
all.equal( coef( sa1iCapLabTen ) * c( 1, 10, 1, 1 ), coef( sa1iCapLab ),
   check.attributes = FALSE, tol = 1e-4 )

sa1iCapTenLabTen <- sfa( mlogOutput ~ logCapitalTen + logLabourTen - 1, 
   data = front41Data )
round( sa1iCapTenLabTen$gridAdj, 3 )
all.equal( coef( sa1iCapTenLabTen, which = "ols" ) * c( 10, 10, 1 ), 
   coef( sa1iCapLab, which = "ols" ), check.attributes = FALSE )
all.equal( coef( sa1iCapTenLabTen, which = "grid" ) * c( 1, 10, 1, 1 ), 
   coef( sa1iCapTenLab, which = "grid" ), check.attributes = FALSE, tol = 1e-4 )
all.equal( coef( sa1iCapTenLabTen ) * c( 10, 10, 1, 1 ), coef( sa1iCapLab ),
   check.attributes = FALSE, tol = 1e-4 )


################################################
## endogenous variable (seemingly) NOT logged ##
################################################

## example data included in FRONTIER 4.1 (cross-section data)
## cross-section data, error components frontier
print( summary( Sa1, logDepVar = FALSE ), digits = 1 )
print( summary( Sa1, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( a1, logDepVar = FALSE ), 2 )
round( efficiencies( a1, asInData = TRUE, logDepVar = FALSE ), 2 )
round( efficiencies( a1, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( a1, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## cross-section data, error components frontier, truncNorm
print( summary( a2, logDepVar = FALSE ), digits = 1 )
round( efficiencies( a2, logDepVar = FALSE ), 2 )
round( efficiencies( a2, asInData = TRUE, logDepVar = FALSE ), 2 )

## cross-section data, efficiency effects frontier
print( summary( Saa1, logDepVar = FALSE ), digits = 1 )
print( summary( Saa1, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
printME( efficiencies( aa1, logDepVar = FALSE, margEff = TRUE ) )
printME( efficiencies( aa1, asInData = TRUE, logDepVar = FALSE ) )
printME( efficiencies( aa1, logDepVar = FALSE, minusU = FALSE ) )
printME( efficiencies( aa1, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ) )

## cross-section data, efficiency effects frontier, zIntercept
print( summary( aa2, logDepVar = FALSE ), digits = 1 )
round( efficiencies( aa2, logDepVar = FALSE ), 2 )
round( efficiencies( aa2, asInData = TRUE, logDepVar = FALSE ), 2 )

## cross-section rice data, error components cost frontier
print( summary( Sdd1, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( dd1, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( dd1, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## cross-section rice data, error components cost frontier, truncNorm
print( summary( dd2, logDepVar = FALSE ), digits = 1 )
print( summary( dd2, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( dd2, logDepVar = FALSE ), 2 )
round( efficiencies( dd2, asInData = TRUE, logDepVar = FALSE ), 2 )
round( efficiencies( dd2, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( dd2, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## cross-section rice data, efficiency effects cost frontier
print( summary( Sdd5, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( dd5, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( dd5, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## cross-section rice data, efficiency effects cost frontier, zIntercept
print( summary( dd6, logDepVar = FALSE ), digits = 1 )
print( summary( dd6, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( dd6, logDepVar = FALSE ), 2 )
round( efficiencies( dd6, asInData = TRUE , logDepVar = FALSE ), 2 )
round( efficiencies( dd6, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( dd6, asInData = TRUE , logDepVar = FALSE, minusU = FALSE ), 2 )

## panel data, error components frontier
print( summary( Sb1, logDepVar = FALSE ), digits = 1 )
round( efficiencies( b1, logDepVar = FALSE ), 2 )
round( efficiencies( b1, asInData = TRUE, logDepVar = FALSE ), 2 )

## panel data, error components frontier, truncNorm
print( summary( b2, logDepVar = FALSE ), digits = 1 )
round( efficiencies( b2, logDepVar = FALSE ), 2 )
round( efficiencies( b2, asInData = TRUE, logDepVar = FALSE ), 2 )

## panel data, error components frontier, timeEffect
print( summary( b3, logDepVar = FALSE ), digits = 1 )
round( efficiencies( b3, logDepVar = FALSE ), 2 )
round( efficiencies( b3, asInData = TRUE, logDepVar = FALSE ), 2 )

## panel data, error components frontier, truncNorm, timeEffect
print( summary( b4, logDepVar = FALSE ), digits = 1 )
print( summary( b4, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( b4, logDepVar = FALSE ), 2 )
round( efficiencies( b4, asInData = TRUE, logDepVar = FALSE ), 2 )
round( efficiencies( b4, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( b4, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## panel data, efficiency effects frontier
print( summary( Sb5, logDepVar = FALSE ), digits = 1 )
round( efficiencies( b5, logDepVar = FALSE ), 2 )
round( efficiencies( b5, asInData = TRUE, logDepVar = FALSE ), 2 )

## panel data, efficiency effects frontier, zIntercept
print( summary( b6, logDepVar = FALSE ), digits = 1 )
print( summary( b6, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( b6, logDepVar = FALSE ), 2 )
round( efficiencies( b6, asInData = TRUE, logDepVar = FALSE ), 2 )
round( efficiencies( b6, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( b6, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## panel rice data, error components cost frontier
print( summary( Sd1, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( d1, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( d1, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## panel rice data, error components cost frontier, truncNorm
print( summary( d2, logDepVar = FALSE ), digits = 1 )
print( summary( d2, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( d2, logDepVar = FALSE ), 2 )
round( efficiencies( d2, asInData = TRUE, logDepVar = FALSE ), 2 )
round( efficiencies( d2, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( d2, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## panel rice data, error components cost frontier, timeEffect
print( summary( d3, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( d3, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( d3, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## panel rice data, error components cost frontier, truncNorm, timeEffect
print( summary( d4, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( d4, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( d4, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## panel rice data, efficiency effects cost frontier
print( summary( Sd5, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( d5, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( d5, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## panel rice data, efficiency effects cost frontier, zIntercept
print( summary( d6, logDepVar = FALSE ), digits = 1 )
print( summary( d6, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( d6, logDepVar = FALSE ), 2 )
round( efficiencies( d6, asInData = TRUE, logDepVar = FALSE ), 2 )
round( efficiencies( d6, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( d6, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## unbalanced panel data, error components frontier
print( summary( b1u, logDepVar = FALSE ), digits = 1 )
round( efficiencies( b1u, logDepVar = FALSE ), 2 )
round( efficiencies( b1u, asInData = TRUE, logDepVar = FALSE ), 2 )

## unbalanced panel data, error components frontier, truncNorm
print( summary( b2u, logDepVar = FALSE ), digits = 1 )
round( efficiencies( b2u, logDepVar = FALSE ), 2 )
round( efficiencies( b2u, asInData = TRUE, logDepVar = FALSE ), 2 )

## unbalanced panel data, error components frontier, timeEffect
print( summary( b3u, logDepVar = FALSE ), digits = 1 )
round( efficiencies( b3u, logDepVar = FALSE ), 2 )
round( efficiencies( b3u, asInData = TRUE, logDepVar = FALSE ), 2 )

## unbalanced panel data, error components frontier, truncNorm, timeEffect
print( summary( b4u, logDepVar = FALSE ), digits = 1 )
round( efficiencies( b4u, logDepVar = FALSE ), 2 )
round( efficiencies( b4u, asInData = TRUE, logDepVar = FALSE ), 2 )

## unbalanced panel data, efficiency effects frontier
print( summary( b5u, logDepVar = FALSE ), digits = 1 )
round( efficiencies( b5u, logDepVar = FALSE ), 2 )
round( efficiencies( b5u, asInData = TRUE, logDepVar = FALSE ), 2 )

## unbalanced panel data, efficiency effects frontier, zIntercept
print( summary( b6u, logDepVar = FALSE ), digits = 1 )
round( efficiencies( b6u, logDepVar = FALSE ), 2 )
round( efficiencies( b6u, asInData = TRUE, logDepVar = FALSE ), 2 )

## unbalanced panel rice data, error components cost frontier
print( summary( d1u, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( d1u, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( d1u, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## unbalanced panel rice data, error components cost frontier, truncNorm
print( summary( d2u, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( d2u, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( d2u, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## unbalanced panel rice data, error components cost frontier, timeEffect
print( summary( d3u, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( d3u, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( d3u, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## unbalanced panel rice data, error components cost frontier, truncNorm, timeEffect
print( summary( d4u, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( d4u, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( d4u, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## unbalanced panel rice data, efficiency effects cost frontier
print( summary( d5u, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( d5u, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( d5u, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )

## unbalanced panel rice data, efficiency effects cost frontier, zIntercept
print( summary( d6u, logDepVar = FALSE, effMinusU = FALSE ), digits = 1 )
round( efficiencies( d6u, logDepVar = FALSE, minusU = FALSE ), 2 )
round( efficiencies( d6u, asInData = TRUE, logDepVar = FALSE, minusU = FALSE ), 2 )


##############################################
## estimation with data NOT in a data frame ##
##############################################

## example data included in FRONTIER 4.1 (cross-section data)
y <- front41Data$output
x1 <- front41Data$capital
x2 <- front41Data$labour
z1 <- front41Data$firmNo

## cross-section data, error components frontier
a1a <- sfa( log( y ) ~ log( x1 ) + log( x2 ) )
all.equal( a1a[-42], a1[-42], check.attributes = FALSE, tol = 1e-4 )
nobs( a1a )

## cross-section data, efficiency effects frontier
aa1a <- sfa( log( y ) ~ log( x1 ) + log( x2 ) | z1 - 1 )
all.equal( aa1a[-42], aa1[-42], check.attributes = FALSE, tol = 1e-4 )

## cross-section data, efficiency effects frontier, zIntercept
aa2a <- sfa( log( y ) ~ log( x1 ) + log( x2 ) | z1 )
all.equal( aa2a[-42], aa2[-42], check.attributes = FALSE, tol = 1e-4 )


##############################################
### estimations with 0 or 1 variable only ###
##############################################

## cross-section data, error components frontier
sa10 <- sfa( logOutput ~ 1, data = front41Data )
a10 <- frontier( "logOutput", NULL, data = front41Data )
print( sa10, digits = 1 )
all.equal( sa10[-42], a10[-42], check.attributes = FALSE, tol = 1e-4 )
nobs( sa10 )

sa11 <- sfa( logOutput ~ logLabour, data = front41Data )
a11 <- frontier( "logOutput", "logLabour", data = front41Data )
print( sa11, digits = 1 )
all.equal( sa11[-42], a11[-42], check.attributes = FALSE, tol = 1e-4 )
nobs( sa11 )

## cross-section data, efficiency effects frontier
saa10 <- sfa( logOutput ~ 1 | firmNo - 1, data = front41Data )
aa10 <- frontier( data = front41Data, "logOutput", NULL,
   zNames = "firmNo" )
print( saa10, digits = 1 )
all.equal( saa10[-42], aa10[-42], tol = 1e-4 )
nobs( saa10 )

saa11 <- sfa( logOutput ~ logLabour | firmNo - 1, data = front41Data )
aa11 <- frontier( data = front41Data, "logOutput", "logLabour",
   zNames = "firmNo" )
print( saa11, digits = 1 )
all.equal( saa11[-42], aa11[-42], tol = 1e-4 )
nobs( saa11 )


##############################################
##### evaluating log likelihood values #######
##############################################
options( digits = 9 )

## cross-section data, error components frontier
print( logLik( a1 ), digits = 4 )
print( logLik( a1, newParam = coef( a1 ) ), digits = 4 )
print( logLik( sa1, newParam = coef( sa1 ) ), digits = 4 )
print( logLik( Sa1, newParam = coef( a1 ) ), digits = 4 )

## cross-section data, error components frontier, truncNorm
print( logLik( a2 ), digits = 4 )
print( logLik( a2, newParam = coef( a2 ) ), digits = 4 )
print( logLik( sa2, newParam = coef( sa2 ) ), digits = 4 )

## cross-section data, error components frontier, truncNorm, starting values
print( logLik( a5 ), digits = 4 )
print( logLik( a5, newParam = coef( a5 ) ), digits = 4 )
print( logLik( sa5, newParam = coef( sa5 ) ), digits = 4 )

## cross-section data, efficiency effects frontier
print( logLik( aa1 ), digits = 4 )
print( logLik( aa1, newParam = coef( aa1 ) ), digits = 4 )
print( logLik( saa1, newParam = coef( saa1 ) ), digits = 4 )
print( logLik( Saa1, newParam = coef( aa1 ) ), digits = 4 )

## cross-section data, efficiency effects frontier, zIntercept
print( logLik( aa2 ), digits = 4 )
print( logLik( aa2, newParam = coef( aa2 ) ), digits = 4 )
print( logLik( saa2, newParam = coef( saa2 ) ), digits = 4 )

## cross-section data, efficiency effects frontier, zIntercept, starting values
print( logLik( aa5 ), digits = 4 )
print( logLik( aa5, newParam = coef( aa5 ) ), digits = 4 )
print( logLik( saa5, newParam = coef( saa5 ) ), digits = 4 )


## data set of rice producers in the Philippines

## cross-section rice data, error components frontier
print( logLik( bb1 ), digits = 4 )
print( logLik( bb1, newParam = coef( bb1 ) ), digits = 4 )
print( logLik( sbb1, newParam = coef( sbb1 ) ), digits = 4 )
print( logLik( Sbb1, newParam = coef( bb1 ) ), digits = 4 )

## cross-section rice data, error components frontier, truncNorm
print( logLik( bb2 ), digits = 4 )
print( logLik( bb2, newParam = coef( bb2 ) ), digits = 4 )
print( logLik( sbb2, newParam = coef( sbb2 ) ), digits = 4 )

## cross-section rice data, efficiency effects frontier
print( logLik( bb5 ), digits = 4 )
print( logLik( bb5, newParam = coef( bb5 ) ), digits = 4 )
print( logLik( sbb5, newParam = coef( sbb5 ) ), digits = 4 )
print( logLik( Sbb5, newParam = coef( bb5 ) ), digits = 4 )

## cross-section rice data, efficiency effects frontier, zIntercept
print( logLik( bb6 ), digits = 4 )
print( logLik( bb6, newParam = coef( bb6 ) ), digits = 4 )
print( logLik( sbb6, newParam = coef( sbb6 ) ), digits = 4 )

## cross-section rice data, error components frontier, truncNorm, starting values
print( logLik( bb7 ), digits = 4 )
print( logLik( bb7, newParam = coef( bb7 ) ), digits = 4 )
print( logLik( sbb7, newParam = coef( sbb7 ) ), digits = 4 )

## cross-section rice data, efficiency effects frontier, zIntercept, starting values
print( logLik( bb8 ), digits = 4 )
print( logLik( bb8, newParam = coef( bb8 ) ), digits = 4 )
print( logLik( sbb8, newParam = coef( sbb8 ) ), digits = 4 )


## Cost Frontier (with land as quasi-fixed input)
## cross-section rice data, error components cost frontier
print( logLik( dd1 ), digits = 4 )
print( logLik( dd1, newParam = coef( dd1 ) ), digits = 4 )
print( logLik( sdd1, newParam = coef( sdd1 ) ), digits = 4 )
print( logLik( Sdd1, newParam = coef( dd1 ) ), digits = 4 )

## cross-section rice data, error components cost frontier, truncNorm
print( logLik( dd2 ), digits = 4 )
print( logLik( dd2, newParam = coef( dd2 ) ), digits = 4 )
print( logLik( sdd2, newParam = coef( sdd2 ) ), digits = 4 )

## cross-section rice data, efficiency effects cost frontier
print( logLik( dd5 ), digits = 4 )
print( logLik( dd5, newParam = coef( dd5 ) ), digits = 4 )
print( logLik( sdd5, newParam = coef( sdd5 ) ), digits = 4 )
print( logLik( Sdd5, newParam = coef( dd5 ) ), digits = 4 )

## cross-section rice data, efficiency effects cost frontier, zIntercept
print( logLik( dd6 ), digits = 4 )
print( logLik( dd6, newParam = coef( dd6 ) ), digits = 4 )
print( logLik( sdd6, newParam = coef( sdd6 ) ), digits = 4 )


## panel data

## panel data, error components frontier
print( logLik( b1 ), digits = 4 )
print( logLik( b1, newParam = coef( b1 ) ), digits = 4 )
print( logLik( sb1, newParam = coef( sb1 ) ), digits = 4 )
print( logLik( Sb1, newParam = coef( b1 ) ), digits = 4 )

## panel data, error components frontier, truncNorm
print( logLik( b2 ), digits = 4 )
print( logLik( b2, newParam = coef( b2 ) ), digits = 4 )
print( logLik( sb2, newParam = coef( sb2 ) ), digits = 4 )

## panel data, error components frontier, timeEffect
print( logLik( b3 ), digits = 4 )
print( logLik( b3, newParam = coef( b3 ) ), digits = 4 )
print( logLik( sb3, newParam = coef( sb3 ) ), digits = 4 )

## panel data, error components frontier, truncNorm, timeEffect
print( logLik( b4 ), digits = 4 )
print( logLik( b4, newParam = coef( b4 ) ), digits = 4 )
print( logLik( sb4, newParam = coef( sb4 ) ), digits = 4 )

## panel data, efficiency effects frontier
print( logLik( b5 ), digits = 4 )
print( logLik( b5, newParam = coef( b5 ) ), digits = 4 )
print( logLik( sb5, newParam = coef( sb5 ) ), digits = 4 )
print( logLik( Sb5, newParam = coef( b5 ) ), digits = 4 )

## panel data, efficiency effects frontier, zIntercept
print( logLik( b6 ), digits = 4 )
print( logLik( b6, newParam = coef( b6 ) ), digits = 4 )
print( logLik( sb6, newParam = coef( sb6 ) ), digits = 4 )

## panel data, error components frontier, truncNorm, timeEffect, starting values
print( logLik( b7 ), digits = 4 )
print( logLik( b7, newParam = coef( b7 ) ), digits = 4 )
print( logLik( sb7, newParam = coef( sb7 ) ), digits = 4 )

## panel data, efficiency effects frontier, zIntercept, starting values
print( logLik( b8 ), digits = 4 )
print( logLik( b8, newParam = coef( b8 ) ), digits = 4 )
print( logLik( sb8, newParam = coef( sb8 ) ), digits = 4 )


## Cost Frontier (with land as quasi-fixed input)
## panel rice data, error components cost frontier
print( logLik( d1 ), digits = 4 )
print( logLik( d1, newParam = coef( d1 ) ), digits = 4 )
print( logLik( sd1, newParam = coef( sd1 ) ), digits = 4 )
print( logLik( Sd1, newParam = coef( d1 ) ), digits = 4 )

## panel rice data, error components cost frontier, truncNorm
print( logLik( d2 ), digits = 4 )
print( logLik( d2, newParam = coef( d2 ) ), digits = 4 )
print( logLik( sd2, newParam = coef( sd2 ) ), digits = 4 )

## panel rice data, error components cost frontier, timeEffect
print( logLik( d3 ), digits = 4 )
print( logLik( d3, newParam = coef( d3 ) ), digits = 4 )
print( logLik( sd3, newParam = coef( sd3 ) ), digits = 4 )

## panel rice data, error components cost frontier, truncNorm, timeEffect
print( logLik( d4 ), digits = 4 )
print( logLik( d4, newParam = coef( d4 ) ), digits = 4 )
print( logLik( sd4, newParam = coef( sd4 ) ), digits = 4 )

## panel rice data, efficiency effects cost frontier
print( logLik( d5 ), digits = 4 )
print( logLik( d5, newParam = coef( d5 ) ), digits = 4 )
print( logLik( sd5, newParam = coef( sd5 ) ), digits = 4 )
print( logLik( Sd5, newParam = coef( d5 ) ), digits = 4 )

## panel rice data, efficiency effects cost frontier, zIntercept
print( logLik( d6 ), digits = 4 )
print( logLik( d6, newParam = coef( d6 ) ), digits = 4 )
print( logLik( sd6, newParam = coef( sd6 ) ), digits = 4 )


## translog frontiers
## cross-section data, error components frontier, translog
print( logLik( translog ), digits = 4 )
print( logLik( translog, newParam = coef( translog ) ), digits = 4 )

## cross-section data, error components frontier, translog, shifter
print( logLik( translogShift ), digits = 4 )
print( logLik( translogShift, newParam = coef( translogShift ) ), digits = 4 )

## cross-section data, efficiency effects frontier, translog
print( logLik( translogZvar ), digits = 4 )
print( logLik( translogZvar, newParam = coef( translogZvar ) ), digits = 4 )


##############################################
########   likelihood ratio tests   ##########
##############################################

## cross-section data, error components frontier
lrtest( a2, a1, a5 )
lrtest( a1, a2, a5 )

## cross-section data, efficiency effects frontier
lrtest( aa2, aa1, aa9, aa2, aa5 )
lrtest( aa9, aa1, aa2, aa5 )

## cross-section data, ECM + EEF
try( lrtest( a2, a1, aa1 ) )
try( lrtest( aa2, a1, aa1 ) )


## data set of rice producers in the Philippines
## cross-section rice data, error components frontier
lrtest( bb2, bb1, bb7 )

## cross-section rice data, efficiency effects frontier
lrtest( bb6, bb5, bb9, bb6, bb8 )


## Cost Frontier (with land as quasi-fixed input)
## cross-section rice data, error components cost frontier
lrtest( dd1, dd2 )

## cross-section rice data, efficiency effects frontier
lrtest( dd6, dd5, dd9, dd6 )


## panel data
## panel data, error components frontier
lrtest( b4, b3, b1, b7 )
lrtest( b4, b2, b1 )
lrtest( b4, b1 )

## panel data, efficiency effects frontier
lrtest( b6, b5, b9, b6, b8 )


## Cost Frontier (with land as quasi-fixed input)
## panel rice data, error components cost frontier
lrtest( d4, d3, d1 )
lrtest( d4, d2, d1 )
lrtest( d4, d1 )

## panel rice data, efficiency effects cost frontier
lrtest( d6, d5, d9, d6 )


## translog
lrtest( translogShift, translog )
