library( "frontier" )
library( "plm" )
options( digits = 5 )

## example data included in FRONTIER 4.1 (cross-section data)
data( front41Data )
front41Data$firmNo <- c( 1:nrow( front41Data ) )

## non-existing variable
try( sfa( log( output ) ~ log( capital7 ) + log( labour ),
   data = front41Data ) )

## nParamTotal > nObs
try( sfa( log( output ) ~ log( capital ) + log( labour ),
   data = front41Data[ 1:4, ] ) )

## nParamTotal >> nObs
try( sfa( log( output ) ~ log( capital ) + log( labour ),
   data = front41Data[ 1:2, ] ) )

## nParamTotal > number of valid observations
try( sfa( log( output ) ~ log( capital ) + log( labour ) + log( firmNo - 56 ),
   data = front41Data ) )

## the dependent variable has only infinite values
try( sfa( log( 0 * output ) ~ log( capital ) + log( labour ),
   data = front41Data ) )

## the dependent variable has only NA values
try( sfa( log( -output ) ~ log( capital ) + log( labour ),
   data = front41Data ) )

## one of the regressors has only infinite values
try( sfa( log( output ) ~ log( 0 * capital ) + log( labour ),
   data = front41Data ) )

## one of the regressors has only NA values
try( sfa( log( output ) ~ log( capital ) + log( -labour ),
   data = front41Data ) )

## one of the regressors of the inefficiency term has only infinite values
try( sfa( log( output ) ~ log( capital ) + log( labour ) | log( 0 * firmNo ),
   data = front41Data ) )

## one of the regressors of the inefficiency term has only NA values
try( sfa( log( output ) ~ log( capital ) + log( labour ) | log(-firmNo ),
   data = front41Data ) )

## no convergence
a1 <- sfa( log( output ) ~ log( capital ) + log( labour ),
   data = front41Data, maxit = 2 )
print( summary( a1 ), digits = 2 )

## no convergence, L(MLE) < L(OLS)
a2 <- sfa( log( output ) ~ log( capital ) + log( labour ),
   data = front41Data, maxit = 2, start = c( 1, 0, 0, 1, 0.5 ) )
print( summary( a2 ), digits = 2 )

## no convergence, L(MLE) < L(OLS), wrong skewness
a3 <- sfa( log( output ) ~ log( capital ) + log( labour ),
   data = front41Data, maxit = 2, ineffDecrease = FALSE )
print( summary( a3, effMinusU = FALSE ), digits = 1 )

## L(MLE) < L(OLS)
a4 <- sfa( log( output ) ~ log( capital ) + log( labour ),
   data = front41Data, start = c( 1, 0, 0, 1, 0.999995 ) )
print( summary( a4 ), digits = 1 )

## too many starting values 
try( sfa( log( output) ~ log( capital ) + log( labour ), data = front41Data,
   truncNorm = TRUE, startVal = c( 0.5, 0.3, 0.5, 0.5, 0.9, -1, 0.3 ) ) )

## too few starting values 
try( sfa( log( output) ~ log( capital ) + log( labour ), data = front41Data,
   truncNorm = TRUE, startVal = c( 0.5, 0.3, 0.5, 0.5, 0.9 ) ) )

## load data abour rice production in the Phillipines
data( "riceProdPhil")

## nobs > nn * nt 
rd <- riceProdPhil
rd <- rbind( rd, rd[ 11, ] )
rd <- plm.data( rd, c( "FMERCODE", "YEARDUM" ) )
try( sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ), data = rd ) )

## non-positive firm number (works now)
rd <- riceProdPhil
rd$FMERCODE <- rd$FMERCODE - 2
rd <- plm.data( rd, c( "FMERCODE", "YEARDUM" ) )
b1 <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ), data = rd )
print( summary( b1 ), digits = 1 )
round( efficiencies( b1 ), 2 )

## firm number > number of firms (works now)
rd <- riceProdPhil
rd$FMERCODE[ rd$FMERCODE == 9 ] <- 47
rd <- plm.data( rd, c( "FMERCODE", "YEARDUM" ) )
b2 <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ), data = rd )
print( summary( b2 ), digits = 1 )
round( efficiencies( b2 ), 2 )
# now with NA
rd <- riceProdPhil
rd$PROD[ rd$FMERCODE == 22 ] <- NA
rd <- plm.data( rd, c( "FMERCODE", "YEARDUM" ) )
b2b <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ), data = rd )
print( summary( b2b ), digits = 1 )
round( efficiencies( b2b ), 2 )

## non-positive period number (works now)
rd <- riceProdPhil
rd$YEARDUM <- rd$YEARDUM - 2
rd <- plm.data( rd, c( "FMERCODE", "YEARDUM" ) )
b3 <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ), data = rd )
print( summary( b3 ), digits = 1 )

## period number > number of periods (works now)
rd <- riceProdPhil
rd$YEARDUM[ rd$YEARDUM == 4 ] <- 10
rd <- plm.data( rd, c( "FMERCODE", "YEARDUM" ) )
b4 <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ), data = rd )
print( summary( b4 ), digits = 1 )
# now with NA
rd <- riceProdPhil
rd$AREA[ rd$YEARDUM == 4 ] <- NA
rd <- plm.data( rd, c( "FMERCODE", "YEARDUM" ) )
b4b <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ), data = rd )
print( summary( b4b ), digits = 1 )

## firm without valid observations (works now)
rd <- riceProdPhil
rd$PROD[ rd$FMERCODE == 12 ] <- NA
rd <- plm.data( rd, c( "FMERCODE", "YEARDUM" ) )
b5 <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ), data = rd )
print( summary( b5 ), digits = 1 )
round( efficiencies( b5 ), 2 )

