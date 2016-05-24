library( frontier )
options( digits = 5 )

## example data included in FRONTIER 4.1 (cross-section data)
data( front41Data )
data( riceProdPhil )
riceProdPhil$cost <- riceProdPhil$LABOR * riceProdPhil$LABORP +
   riceProdPhil$NPK * riceProdPhil$NPKP

###### left-skewed residuals but ineffDecrease = FALSE
## front41Data
a1 <- sfa( log( output ) ~ log( capital ) + log( labour ),
   ineffDecrease = FALSE, data = front41Data )
print( summary( a1, effMinusU = FALSE ), digits = 1 )
lrtest( a1 )

## front41Data, truncNorm
a2 <- sfa( log( output ) ~ log( capital ) + log( labour ),
   ineffDecrease = FALSE, truncNorm = TRUE, data = front41Data )
print( summary( a2, effMinusU = FALSE ), digits = 1 )
lrtest( a2 )

## riceProdPhil
b1 <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ),
   ineffDecrease = FALSE, data = riceProdPhil )
print( summary( b1, effMinusU = FALSE ), digits = 1 )
lrtest( b1 )

## riceProdPhil, truncNorm
b2 <- sfa( log( PROD ) ~ log( AREA ) + log( LABOR ) + log( NPK ),
   ineffDecrease = FALSE, truncNorm = TRUE, data = riceProdPhil )
print( summary( b2, effMinusU = FALSE ), digits = 1 )
lrtest( b2 )


###### right-skewed residuals but ineffDecrease = TRUE
## riceProdPhil
d1 <- sfa( log( cost ) ~ log( PROD ) + log( AREA ) + log( LABORP ) +
   log( NPKP ), data = riceProdPhil )
print( summary( d1 ), digits = 1 )
lrtest( d1 )

## riceProdPhil, truncNorm
d2 <- sfa( log( cost ) ~ log( PROD ) + log( AREA ) + log( LABORP ) +
   log( NPKP ), truncNorm = TRUE, data = riceProdPhil )
print( summary( d2 ), digits = 1 )
lrtest( d2 )
