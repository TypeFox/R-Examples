library( micEcon )
options( digits = 3 )

## preparing data
data( germanFarms )
# quantity of variable inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# quantity of crop outputs
germanFarms$qCrop <- germanFarms$vCrop / germanFarms$pOutput
# quantity of animal outputs
germanFarms$qAnimal <- germanFarms$vAnimal / germanFarms$pOutput
# a time trend to account for technical progress:
germanFarms$time <- c(1:20)
# generate (artificial) prices
germanFarms$pLand <- 200 + 15 * germanFarms$time

# estimate a translog ray production function
estResultRay <- translogRayEst( yNames = c( "qCrop", "qAnimal" ),
   xNames = c( "qLabor", "land", "qVarInput" ),
   data = germanFarms )
print( estResultRay )
summary( estResultRay )
print.default( estResultRay )

# different order of outputs
estResultRay2 <- translogRayEst( yNames = c( "qAnimal", "qCrop" ),
   xNames = c( "qLabor", "land", "qVarInput" ),
   data = germanFarms )
print( estResultRay2 )
summary( estResultRay2 )
all.equal( abs( coef( estResultRay2 )[ 6:15 ] ),
   abs( coef( estResultRay )[ 6:15 ] ) )

# different order of inputs
estResultRay3 <- translogRayEst( yNames = c( "qCrop", "qAnimal" ),
   xNames = c( "qVarInput", "qLabor", "land" ),
   data = germanFarms )
print( estResultRay3 )
summary( estResultRay3 )
all.equal( coef( estResultRay ), coef( estResultRay3 )[
   c( 1, 3, 4, 2, 5, 10, 11, 7, 12, 13, 8, 14, 6, 9, 15 ) ],
   check.attributes = FALSE )


## testing translogRayDeriv
tlRayDeriv <- translogRayDeriv( yNames = c( "qCrop", "qAnimal" ),
   xNames = c( "qLabor", "land", "qVarInput" ),
   data = germanFarms, coef = coef( estResultRay ) )
print( tlRayDeriv )

tlRayDeriv2 <- translogRayDeriv( yNames = c( "qAnimal", "qCrop" ),
   xNames = c( "qLabor", "land", "qVarInput" ),
   data = germanFarms, coef = coef( estResultRay2 ) )
all.equal( tlRayDeriv, tlRayDeriv2[ , c( 1:3, 5, 4 ) ] )

tlRayDeriv3 <- translogRayDeriv( yNames = c( "qCrop", "qAnimal" ),
   xNames = c( "qVarInput", "qLabor", "land" ),
   data = germanFarms, coef = coef( estResultRay3 ) )
all.equal( tlRayDeriv, tlRayDeriv3[ , c( 2, 3, 1, 4, 5 ) ] )


## testing translogProdFuncMargCost with a ray function
# compute the marginal costs of producing the output
margCostRay <- translogProdFuncMargCost( yNames = c( "qCrop", "qAnimal" ),
   xNames = c( "qLabor", "land", "qVarInput" ),
   wNames = c( "pLabor", "pLand", "pVarInput" ),
   data = germanFarms, coef = coef( estResultRay ) )
print( margCostRay )

margCostRay2 <- translogProdFuncMargCost( yNames = c( "qAnimal", "qCrop" ),
   xNames = c( "qLabor", "land", "qVarInput" ),
   wNames = c( "pLabor", "pLand", "pVarInput" ),
   data = germanFarms, coef = coef( estResultRay2 ) )
all.equal( margCostRay, margCostRay2[ , c( 2:1 ) ] )

margCostRay3 <- translogProdFuncMargCost( yNames = c( "qCrop", "qAnimal" ),
   xNames = c( "qVarInput", "qLabor", "land" ),
   wNames = c( "pVarInput", "pLabor", "pLand" ),
   data = germanFarms, coef = coef( estResultRay3 ) )
all.equal( margCostRay, margCostRay3 )


######################################
##   using data set appleProdFr86   ##
######################################
data( "appleProdFr86" )
# quantity of the capital input
appleProdFr86$qCap <- with( appleProdFr86, vCap / pCap )
# quantity of the labour input
appleProdFr86$qLab <- with( appleProdFr86, vLab / pLab )
# quantity of the materials input
appleProdFr86$qMat <- with( appleProdFr86, vMat / pMat )

# estimate a translog ray production function
estApple <- translogRayEst( yNames = c( "qApples", "qOtherOut" ),
   xNames = c( "qCap", "qLab", "qMat" ), data = appleProdFr86 )
print( estApple )
summary( estApple )
print.default( estApple )

# different order of outputs
estApple2 <- translogRayEst( yNames = c( "qOtherOut", "qApples" ),
   xNames = c( "qCap", "qLab", "qMat" ), data = appleProdFr86 )
print( estApple2 )
summary( estApple2 )
all.equal( abs( coef( estApple2 )[ 6:15 ] ),
   abs( coef( estApple )[ 6:15 ] ) )

# different order of inputs
estApple3 <- translogRayEst( yNames = c( "qApples", "qOtherOut" ),
   xNames = c( "qMat", "qCap", "qLab" ), data = appleProdFr86 )
print( estApple3 )
summary( estApple3 )
all.equal( coef( estApple ), coef( estApple3 )[
   c( 1, 3, 4, 2, 5, 10, 11, 7, 12, 13, 8, 14, 6, 9, 15 ) ],
   check.attributes = FALSE )


## testing translogRayDeriv
derivApple <- translogRayDeriv( yNames = c( "qApples", "qOtherOut" ),
   xNames = c( "qCap", "qLab", "qMat" ), data = appleProdFr86,
   coef = coef( estApple ) )
print( derivApple )

derivApple2 <- translogRayDeriv( yNames = c( "qOtherOut", "qApples" ),
   xNames = c( "qCap", "qLab", "qMat" ), data = appleProdFr86,
   coef = coef( estApple2 ) )
all.equal( derivApple, derivApple2[ , c( 1:3, 5, 4 ) ] )

derivApple3 <- translogRayDeriv( yNames = c( "qApples", "qOtherOut" ),
   xNames = c( "qMat", "qCap", "qLab" ), data = appleProdFr86, 
   coef = coef( estApple3 ) )
all.equal( derivApple, derivApple3[ , c( 2, 3, 1, 4, 5 ) ] )


## testing translogProdFuncMargCost with a ray function
# compute the marginal costs of producing the output
margCostApple <- translogProdFuncMargCost( 
   yNames = c( "qApples", "qOtherOut" ),
   xNames = c( "qCap", "qLab", "qMat" ),
   wNames = c( "pCap", "pLab", "pMat" ),
   data = appleProdFr86, coef = coef( estApple ) )
print( margCostApple )

margCostApple2 <- translogProdFuncMargCost( 
   yNames = c( "qOtherOut", "qApples" ),
   xNames = c( "qCap", "qLab", "qMat" ),
   wNames = c( "pCap", "pLab", "pMat" ),
   data = appleProdFr86, coef = coef( estApple2 ) )
all.equal( margCostApple, margCostApple2[ , c( 2:1 ) ] )

margCostApple3 <- translogProdFuncMargCost( 
   yNames = c( "qApples", "qOtherOut" ),
   xNames = c( "qMat", "qCap", "qLab" ),
   wNames = c( "pMat", "pCap", "pLab" ),
   data = appleProdFr86, coef = coef( estApple3 ) )
all.equal( margCostApple, margCostApple3 )
