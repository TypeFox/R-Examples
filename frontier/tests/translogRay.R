library( frontier )
options( digits = 5 )

printAll <- function( x ) {
   for( n in names( x ) ) {
      cat( "$", n, "\n", sep = "" )
      if( n %in% c( "olsParam", "gridParam", "mleParam", "olsStdEr", "mleCov" ) ) {
         print( round( x[[ n ]], 2 ) )
      } else if( n %in% c( "fitted", "resid", "olsResid" ) ) {
         print( round( x[[ n ]], 3 ) )
      } else {
         print( x[[ n ]] )
      }
      cat( "\n" )
   }
   cat( "class\n" )
   print( class( x ) )
}

## preparing data
data( germanFarms )
# quantity of crop outputs
germanFarms$qCrop <- germanFarms$vCrop / germanFarms$pOutput
# quantity of animal outputs
germanFarms$qAnimal <- germanFarms$vAnimal / germanFarms$pOutput
# quantity of variable inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# a time trend to account for technical progress:
germanFarms$time <- c(1:20)

# estimate a translog ray production function
estResultRay <- frontierTranslogRay( yNames = c( "qCrop", "qAnimal" ),
   xNames = c( "qLabor", "land", "qVarInput" ),
   data = germanFarms )
print( estResultRay )
print( summary( estResultRay ), digits = 1 )
lrtest( estResultRay )
round( efficiencies( estResultRay ), 3 )
round( efficiencies( estResultRay, asInData = TRUE ), 3 )
round( residuals( estResultRay ), 3 )
round( residuals( estResultRay, asInData = TRUE ), 3 )
printAll( estResultRay )
