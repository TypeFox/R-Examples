library( micEcon )
options( digits = 3 )

data( germanFarms )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# value of labor input
germanFarms$vLabor <- germanFarms$pLabor + germanFarms$qLabor
# total variable cost
germanFarms$cost <- germanFarms$vLabor + germanFarms$vVarInput
# a time trend to account for technical progress:
germanFarms$time <- c(1:20)

# non-hom in prices, without land and trend
estResult <- translogCostEst( cName = "cost", yName = "qOutput", 
   pNames = c( "pLabor", "pVarInput" ),
   data = germanFarms, homPrice = FALSE )
print( estResult )
summary( estResult$est )

# non-hom in prices, without land, with trend
estResultTrend <- translogCostEst( cName = "cost", yName = "qOutput", 
   pNames = c( "pLabor", "pVarInput" ),
   shifterNames = "time", data = germanFarms, homPrice = FALSE )
print( estResultTrend )
summary( estResultTrend$est )

# non-hom in prices, with land, without trend
estResultLand <- translogCostEst( cName = "cost", yName = "qOutput", 
   pNames = c( "pLabor", "pVarInput" ), fNames = "land",
   data = germanFarms, homPrice = FALSE )
print( estResultLand )
summary( estResultLand$est )

# non-hom in prices, with land + trend
estResultLandTrend <- translogCostEst( cName = "cost", yName = "qOutput", 
   pNames = c( "pLabor", "pVarInput" ), fNames = "land",
   shifterNames = "time", data = germanFarms, homPrice = FALSE )
print( estResultLandTrend )
summary( estResultLandTrend$est )
