library( micEcon )
options( digits = 4 )

## function for printIndexing indices
printIndices <- function( what, ... ) {
   for( i in c( "Laspeyres", "Paasche", "Fisher" ) ) {
      cat( "\n", i, "\n" )
      if( what == "p" ) {
         index <- priceIndex( ..., method = i, weights = TRUE )
      } else {
         index <- quantityIndex( ..., method = i, weights = TRUE )
      }
      print( index )
      testRowSums <- rowSums( attributes( index )$weights[ !is.na( index ), ] )
      names( testRowSums ) <- NULL
      if( all.equal( testRowSums,
         rep( 1, sum( !is.na( index ) ) ) ) != TRUE ) {
         cat( "\nrowSums are not equal to one!!!\n\n" )
      }
   }
}

## Missong03E7.7
data( Missong03E7.7 )
## price indices for Missong03E7.7
printIndices( "p",  c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong03E7.7 )

## quantity indices for Missong03E7.7
printIndices( "q",  c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong03E7.7 )


## Bleymueller79E25.1
data( Bleymueller79E25.1 )
## price indices for Bleymueller79E25.1
printIndices( "p",  c( "p.A", "p.B", "p.C", "p.D" ),
   c( "q.A", "q.B", "q.C", "q.D" ),  1, Bleymueller79E25.1 )

## quantity indices for Bleymueller79E25.1
printIndices( "q",  c( "p.A", "p.B", "p.C", "p.D" ),
   c("q.A", "q.B", "q.C", "q.D" ), 1, Bleymueller79E25.1 )


## Electricity (Christensen & Greene 1976)
data( "Electricity", package = "Ecdat" )
Electricity <- Electricity[ 1:35, ]
## preparing electricity data
pNames <- c( "pl", "pk", "pf" )
qNames <- c( "ql", "qk", "qf" )
sNames <- c( "sl", "sk", "sf" )
for( i in 1:3 ) {
   Electricity[[ qNames[ i ] ]] <- Electricity[[ sNames[ i ] ]] *
      Electricity[[ "cost" ]] / Electricity[[ pNames[ i ] ]]
}
allObs <- !is.na( Electricity$pl )

## price indices for Electricity data
printIndices( "p",  pNames, qNames, 1, Electricity )

## quantity indices for Electricity data
printIndices( "q",  pNames, qNames, 1, Electricity )

## price indices for Electricity data and base=mean
printIndices( "p",  pNames, qNames, allObs, Electricity )

## quantity indices for Electricity data and base=mean
printIndices( "q",  pNames, qNames, allObs, Electricity )


## Electricity data with some NA prices
## manipulating data of Electricity data
ElectricityNaP <- Electricity
for( i in 1:3 ) {
   ElectricityNaP[[ pNames[ i ] ]][ c( 2, i * 4, i * 8 ) ] <- NA
}

## price indices for Electricity data with some NA prices
printIndices( "p",  pNames, qNames, 1, ElectricityNaP )

## quantity indices for Electricity data with some NA prices
printIndices( "q",  pNames, qNames, 1, ElectricityNaP )

## price indices for Electricity data with some NA prices and na.rm=TRUE
printIndices( "p",  pNames, qNames, 1, ElectricityNaP, na.rm = TRUE )

## quantity indices for Electricity data with some NA prices and na.rm=TRUE
printIndices( "q",  pNames, qNames, 1, ElectricityNaP, na.rm = TRUE )

## price indices for Electricity data with some NA prices and base=mean
printIndices( "p",  pNames, qNames, 16, ElectricityNaP )

## quantity indices for Electricity data with some NA prices and base=mean
printIndices( "q",  pNames, qNames, allObs, ElectricityNaP )

## price indices for Electricity data with some NA prices and na.rm=TRUE and base=mean
printIndices( "p",  pNames, qNames, allObs, ElectricityNaP, na.rm = TRUE )

## quantity indices for Electricity data with some NA prices and na.rm=TRUE and base=mean
printIndices( "q",  pNames, qNames, allObs, ElectricityNaP, na.rm = TRUE )


## Electricity data with some NA prices and quantities
## manipulating data of Electricity data
ElectricityNaPQ <- ElectricityNaP
for( i in 1:3 ) {
   ElectricityNaPQ[[ qNames[ i ] ]][ c( 2, ( i + 1 ) * 4, i * 8 ) ] <- NA
}

## price indices for Electricity data with some NA prices and quantities
printIndices( "p",  pNames, qNames, 1, ElectricityNaPQ )

## quantity indices for Electricity data with some NA prices and quantities
printIndices( "q",  pNames, qNames, 1, ElectricityNaPQ )

## price indices for Electricity data with some NA prices and quantities and na.rm=TRUE
printIndices( "p",  pNames, qNames, 1, ElectricityNaPQ, na.rm = TRUE )

## quantity indices for Electricity data with some NA prices and quantities and na.rm=TRUE
printIndices( "q",  pNames, qNames, 1, ElectricityNaPQ, na.rm = TRUE )

## price indices for Electricity data with some NA prices and quantities and base=mean
printIndices( "p",  pNames, qNames, 16, ElectricityNaPQ )

## quantity indices for Electricity data with some NA prices and quantities and base=mean
printIndices( "q",  pNames, qNames, 16, ElectricityNaPQ )

## price indices for Electricity data with some NA prices and quantities and na.rm=TRUE and base=mean
printIndices( "p",  pNames, qNames, allObs, ElectricityNaPQ, na.rm = TRUE )

## quantity indices for Electricity data with some NA prices and quantities and na.rm=TRUE and base=mean
printIndices( "q",  pNames, qNames, allObs, ElectricityNaPQ, na.rm = TRUE )


## Electricity data with some NA prices, where quantities are partly zero
## manipulating data of Electricity data
ElectricityNaP0Q <- Electricity
for( i in 1:3 ) {
   ElectricityNaP0Q[[ pNames[ i ] ]][ c( 2, i * 4, i * 8 ) ] <- NA
   ElectricityNaP0Q[[ qNames[ i ] ]][ c( 2, i * 8 ) ] <- 0
}

## price indices, base = 1
printIndices( "p",  pNames, qNames, 1, ElectricityNaP0Q )

## quantity indices, base = 1
printIndices( "q",  pNames, qNames, 1, ElectricityNaP0Q )

## price indices, base = mean, na.rm = TRUE
printIndices( "p",  pNames, qNames, allObs, ElectricityNaP0Q, na.rm = TRUE )

## quantity indices, base = mean, na.rm = TRUE
printIndices( "q",  pNames, qNames, allObs, ElectricityNaP0Q, na.rm = TRUE )
