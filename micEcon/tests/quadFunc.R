library( micEcon )
library( plm )
options( digits = 3 )
printQuadFuncEst <- function( x ) {
   for( n in names( x ) ) {
      cat( "$", n, "\n", sep ="" )
      if( n == "model.matrix" && is.numeric( x[[ n ]] ) ) {
         for( i in 1:ncol( x[[ n ]] ) ) {
            size <- floor( log( x[[ n ]][ , i ], base = 10 ) ) + 1
            size[ size < -5 ] <- -5
            fac <- 10^( size - 5 )
            x[[ n ]][ , i ] <- fac * round( x[[ n ]][ , i ] / fac ) 
         }
      }
      print( x[[ n ]] )
      cat( "\n" )
   }
   cat( "attr(,\"class\")\n" )
   print( class( x ) )   
}

## preparing data
data( germanFarms )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of variable inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# a time trend to account for technical progress:
germanFarms$time <- c(1:20)
# mean values of all variables
germanFarmsMeans <- colMeans( germanFarms[ , -1 ] )

## estimate a quadratic production function
estResult <- quadFuncEst( "qOutput",
   c( "qLabor", "land", "qVarInput", "time" ), germanFarms )
coef( estResult )
print( estResult )

## compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ) )
all.equal( fitted, estResult$fitted )
# only at mean values
fittedMean <- quadFuncCalc(
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   data = germanFarmsMeans, coef = coef( estResult ) )
print( fittedMean )

## compute the marginal products of the inputs
margProducts <- quadFuncDeriv(
   c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), vcov( estResult ) )
print( margProducts )
print( attributes( margProducts )$variance )
print( attributes( margProducts )$stdDev )
# only at mean values
margProductsMean <- quadFuncDeriv(
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   data = germanFarmsMeans, coef = coef( estResult ),
   coefCov = vcov( estResult ) )
print( margProductsMean )
print( attributes( margProductsMean )$variance )
print( attributes( margProductsMean )$stdDev )

## estimate a quadratic production function with a shifter
estResultShifter <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = "time", data = germanFarms )
coef( estResultShifter )
print( estResultShifter )
# compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = "time", data = germanFarms, coef( estResultShifter ) )
all.equal( fitted, estResultShifter$fitted )
# compute marginal products = partial derivatives
margProdShifter <- quadFuncDeriv(
   c( "qLabor", "land", "qVarInput" ),
   germanFarms, coef( estResultShifter ), vcov( estResultShifter ) )
print( margProdShifter )
print( attributes( margProdShifter )$variance )
print( attributes( margProdShifter )$stdDev )

## estimate a quadratic production function with 2 shifters
germanFarms$timeSq <- germanFarms$time^2
estResultShifter2 <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "time", "timeSq" ), data = germanFarms )
coef( estResultShifter2 )
print( estResultShifter2 )
# compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "time", "timeSq" ), data = germanFarms, 
   coef( estResultShifter2 ) )
all.equal( fitted, estResultShifter2$fitted )

## estimate a linear functions with quadFuncEst
estResultLinear <- quadFuncEst( yName = "qOutput", xNames = NULL,
   shifterNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms )
coef( estResultLinear )
print( estResultLinear )
estResultLin <- quadFuncEst( yName = "qOutput", 
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms, linear = TRUE )
all.equal( coef( estResultLinear ), coef( estResultLin )[1:5],
   check.attributes = FALSE )
all.equal( vcov( estResultLinear ), vcov( estResultLin )[1:5,1:5],
   check.attributes = FALSE )
coef( estResultLin )
vcov( estResultLin )
print( estResultLin )
# compute fitted values
fitted <- quadFuncCalc( xNames = NULL,
   shifterNames = c( "time", "qLabor", "land", "qVarInput" ), 
   data = germanFarms, coef( estResultLinear ) )
all.equal( fitted, estResultLinear$fitted )
fitted <- quadFuncCalc( xNames = c( "time", "qLabor", "land", "qVarInput" ), 
   data = germanFarms, coef( estResultLin ) )
all.equal( fitted, estResultLin$fitted )
all.equal( estResultLinear$fitted, estResultLin$fitted )
# compute partial derivatives
margProducts <- quadFuncDeriv(
   c( "qLabor", "land", "qVarInput", "time" ),
   data = germanFarms, coef = coef( estResultLin ), 
   coefCov = vcov( estResultLin ) )
sapply( margProducts, sd )
all.equal( margProducts[1,], coef( estResultLin )[2:5],
   check.attributes = FALSE )
all.equal( attributes( margProducts )$variance[1,], diag( vcov( estResultLin ) )[2:5], 
   check.attributes = FALSE )

## estimate a quadratic production function with a logical variable as shifter
germanFarms$reUnif <- germanFarms$time >= 16
estResultShifterLogi <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "reUnif" ), data = germanFarms )
coef( estResultShifterLogi )
print( estResultShifterLogi )
# compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "reUnif" ), data = germanFarms, 
   coef( estResultShifterLogi ) )
all.equal( fitted, estResultShifterLogi$fitted )

## estimate a quadratic production function with a factor as shifter
germanFarms$decade <- as.factor( c( rep( "70s", 5 ), rep( "80s", 10 ), 
   rep( "90s", 5 ) ) )
estResultShifterFac <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "decade" ), data = germanFarms )
coef( estResultShifterFac )
print( estResultShifterFac )
# compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "decade" ), data = germanFarms, 
   coef( estResultShifterFac ) )
all.equal( fitted, estResultShifterFac$fitted )

## estimate a quadratic production function with some shifters are logical
estResultShifterLogi2 <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "time", "reUnif" ), data = germanFarms )
coef( estResultShifterLogi2 )
print( estResultShifterLogi2 )
# compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "time", "reUnif" ), data = germanFarms, 
   coef( estResultShifterLogi2 ) )
all.equal( fitted, estResultShifterLogi2$fitted )

## estimate a quadratic production function with some shifters are factors
estResultShifterFac2 <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "time", "decade" ), data = germanFarms )
coef( estResultShifterFac2 )
print( estResultShifterFac2 )
# compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "time", "decade" ), data = germanFarms, 
   coef( estResultShifterFac2 ) )
all.equal( fitted, estResultShifterFac2$fitted )

## estimate with further argument passed to lm()
estResult2 <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, x = TRUE, y = TRUE )
coef( estResult2 )
print( estResult2 )

## calculating elasticities
estElaFit <- quadFuncEla( 
   xNames = c( "qLabor", "land", "qVarInput", "time" ), 
   data = germanFarms, coef = coef( estResult ) )
all.equal( estElaFit, elas( estResult ) )
print( estElaFit )
estElaObs <- quadFuncEla( 
   xNames = c( "qLabor", "land", "qVarInput", "time" ), 
   data = germanFarms, coef = coef( estResult ),
   yName = "qOutput" )
all.equal( estElaObs, elas( estResult, yObs = TRUE ) )
print( estElaObs )
max( abs( estElaFit - estElaObs ) )
# only at mean values
estElaMeanFit <- elas( estResult, data = germanFarmsMeans )
print( estElaMeanFit )
estElaMeanObs <- elas( estResult, data = germanFarmsMeans, yObs = TRUE )
print( estElaMeanObs )
print( estElaMeanFit - estElaMeanObs )
# with a shifter
estElaShifterFit <- quadFuncEla( 
   xNames = c( "qLabor", "land", "qVarInput" ), 
   data = germanFarms, coef = coef( estResultShifter ),
   shifterNames = "time" )
all.equal( estElaShifterFit, elas( estResultShifter ) )
print( estElaShifterFit )
estElaShifterObs <- quadFuncEla( 
   xNames = c( "qLabor", "land", "qVarInput" ), 
   data = germanFarms, coef = coef( estResultShifter ),
   yName = "qOutput" )
all.equal( estElaShifterObs, elas( estResultShifter, yObs = TRUE ) )
print( estElaShifterObs )
max( abs( estElaShifterFit - estElaShifterObs ) )
estElaShifterObs2 <- quadFuncEla( 
   xNames = c( "qLabor", "land", "qVarInput" ), 
   data = germanFarms, coef = coef( estResultShifter ),
   yName = "qOutput", shifterNames = "time" )
all.equal( estElaShifterObs, estElaShifterObs2 )
# only at mean values
estElaShifterMeanFit <- elas( estResultShifter, data = germanFarmsMeans )
print( estElaShifterMeanFit )
estElaShifterMeanObs <- elas( estResultShifter, data = germanFarmsMeans,
   yObs = TRUE )
print( estElaShifterMeanObs )
print( estElaShifterMeanFit - estElaShifterMeanObs )


################ imposing homogeneity #####################
## linear functions with homogeneity imposed
estResultLinHom <- quadFuncEst( yName = "qOutput", 
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms, linear = TRUE, 
   homWeights = c( qLabor = 0.2, land = 0.5, qVarInput = 0.3 ) )
coef( estResultLinHom )
all.equal( sum( coef( estResultLinHom )[3:4] ), - coef( estResultLinHom )[5],
   check.attributes = FALSE )
vcov( estResultLinHom )
all.equal( rowSums( vcov( estResultLinHom )[ , 3:4 ] ), 
   - vcov( estResultLinHom )[ , 5 ] )
all.equal( colSums( vcov( estResultLinHom )[ 3:4, ] ), 
   - vcov( estResultLinHom )[ 5, ] )
# different order of weights
estResultLinHom2 <- quadFuncEst( yName = "qOutput", 
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms, linear = TRUE, 
   homWeights = c( qVarInput = 0.3, land = 0.5, qLabor = 0.2 ) )
all.equal( coef( estResultLinHom ), coef( estResultLinHom2 ) )
all.equal( vcov( estResultLinHom ), vcov( estResultLinHom2 ) )
# different order of xNames
estResultLinHom3 <- quadFuncEst( yName = "qOutput", 
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   data = germanFarms, linear = TRUE, 
   homWeights = c( qLabor = 0.2, land = 0.5, qVarInput = 0.3 ) )
all.equal( coef( estResultLinHom ), 
   coef( estResultLinHom3 )[ c( 1, 5, 2:4, 6:15 ) ],
   check.attributes = FALSE )
all.equal( vcov( estResultLinHom ), 
   vcov( estResultLinHom3 )[ c( 1, 5, 2:4, 6:15 ), c( 1, 5, 2:4, 6:15 ) ],
   check.attributes = FALSE )
# homogenous in all variables
estResultLinHom4 <- quadFuncEst( yName = "qOutput", 
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms, linear = TRUE, 
   homWeights = c( qLabor = 0.2, land = 0.5, qVarInput = 0.3, time = 0 ) )
coef( estResultLinHom4 )
all.equal( sum( coef( estResultLinHom4 )[2:4] ), 
 - coef( estResultLinHom4 )[5], check.attributes = FALSE )
vcov( estResultLinHom4 )
all.equal( rowSums( vcov( estResultLinHom4 )[ , 2:4 ] ), 
   - vcov( estResultLinHom4 )[ , 5 ] )
all.equal( colSums( vcov( estResultLinHom4 )[ 2:4, ] ), 
   - vcov( estResultLinHom4 )[ 5, ] )

## computing fitted values of linear functions with homogeneity imposed
fitted <- quadFuncCalc(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultLinHom ), data = germanFarms, 
   homWeights = c( qLabor = 0.2, land = 0.5, qVarInput = 0.3 ) )
all.equal( estResultLinHom$fitted, fitted )
fitted <- quadFuncCalc(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultLinHom ), data = germanFarms, 
   homWeights = c( qVarInput = 0.3, land = 0.5, qLabor = 0.2 ) )
all.equal( estResultLinHom$fitted, fitted )
fitted <- quadFuncCalc(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultLinHom2 ), data = germanFarms, 
   homWeights = c( qVarInput = 0.3, land = 0.5, qLabor = 0.2 ) )
all.equal( estResultLinHom2$fitted, fitted )
fitted <- quadFuncCalc(
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   coef = coef( estResultLinHom3 ), data = germanFarms,
   homWeights = c( qLabor = 0.2, land = 0.5, qVarInput = 0.3 ) )
all.equal( estResultLinHom3$fitted, fitted )
fitted <- quadFuncCalc(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultLinHom4 ), data = germanFarms,
   homWeights = c( qLabor = 0.2, land = 0.5, qVarInput = 0.3, time = 0 ) )
all.equal( estResultLinHom4$fitted, fitted )

## derivatives of linear functions with homogeneity imposed
estResultLinHomDeriv <- quadFuncDeriv(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultLinHom ), data = germanFarms, 
   homWeights = c( qLabor = 0.2, land = 0.5, qVarInput = 0.3 ) )
print( estResultLinHomDeriv )
all.equal( estResultLinHomDeriv$qLabor * germanFarms$qLabor +
   estResultLinHomDeriv$land * germanFarms$land,
   - estResultLinHomDeriv$qVarInput * germanFarms$qVarInput )
# different order of weights (different from order used for estimation)
estResultLinHom1Deriv <- quadFuncDeriv(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultLinHom ), data = germanFarms, 
   homWeights = c( qVarInput = 0.3, land = 0.5, qLabor = 0.2 ) )
all.equal( estResultLinHomDeriv, estResultLinHom1Deriv )
# different order of weights (same order as used for estimation)
estResultLinHom2Deriv <- quadFuncDeriv(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultLinHom2 ), data = germanFarms, 
   homWeights = c( qVarInput = 0.3, land = 0.5, qLabor = 0.2 ) )
all.equal( estResultLinHomDeriv, estResultLinHom2Deriv )
# different order of independent variables
estResultLinHom3Deriv <- quadFuncDeriv(
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   coef = coef( estResultLinHom3 ), data = germanFarms,
   homWeights = c( qLabor = 0.2, land = 0.5, qVarInput = 0.3 ) )
all.equal( estResultLinHomDeriv, 
   estResultLinHom3Deriv[ , c( 4, 1, 2, 3 ) ] )
# homogenous in all independent variables
estResultLinHom4Deriv <- quadFuncDeriv(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultLinHom4 ), data = germanFarms,
   homWeights = c( qLabor = 0.2, land = 0.5, qVarInput = 0.3, time = 0 ) )
print( estResultLinHom4Deriv )
all.equal( estResultLinHom4Deriv$qLabor * germanFarms$qLabor +
   estResultLinHom4Deriv$land * germanFarms$land +
   estResultLinHom4Deriv$qVarInput * germanFarms$qVarInput,
   - estResultLinHom4Deriv$time * germanFarms$time )

## quadratic functions with homogeneity imposed
estResultHom <- quadFuncEst( yName = "qOutput", 
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms,
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2 ) )
coef( estResultHom )
all.equal( sum( coef( estResultHom )[3:4] ), - coef( estResultHom )[5],
   check.attributes = FALSE ) # a_2:4
all.equal( sum( coef( estResultHom )[7:8] ), - coef( estResultHom )[9],
   check.attributes = FALSE ) # b_1_2:4
all.equal( sum( coef( estResultHom )[10:11] ), - coef( estResultHom )[12],
   check.attributes = FALSE ) # b_2_2:4
all.equal( sum( coef( estResultHom )[13:14] ), - coef( estResultHom )[11],
   check.attributes = FALSE ) # b_3_2:4
all.equal( sum( coef( estResultHom )[14:15] ), - coef( estResultHom )[12],
   check.attributes = FALSE ) # b_4_2:4
vcov( estResultHom )
all.equal( rowSums( vcov( estResultHom )[ , 3:4 ] ), 
   - vcov( estResultHom )[ , 5 ] ) # a_2:4
all.equal( rowSums( vcov( estResultHom )[ , 7:8 ] ), 
   - vcov( estResultHom )[ , 9 ] ) # b_1_2:4
all.equal( rowSums( vcov( estResultHom )[ , 10:11 ] ), 
   - vcov( estResultHom )[ , 12 ] ) # b_2_2:4
all.equal( rowSums( vcov( estResultHom )[ , 13:14 ] ), 
   - vcov( estResultHom )[ , 11 ] ) # b_3_2:4
all.equal( rowSums( vcov( estResultHom )[ , 14:15 ] ), 
   - vcov( estResultHom )[ , 12 ] ) # b_4_2:4
all.equal( colSums( vcov( estResultHom )[ 3:4, ] ), 
   - vcov( estResultHom )[ 5, ] ) # a_2:4
all.equal( colSums( vcov( estResultHom )[ 7:8, ] ), 
   - vcov( estResultHom )[ 9, ] ) # b_1_2:4
all.equal( colSums( vcov( estResultHom )[ 10:11, ] ), 
   - vcov( estResultHom )[ 12, ] ) # b_2_2:4
all.equal( colSums( vcov( estResultHom )[ 13:14, ] ), 
   - vcov( estResultHom )[ 11, ] ) # b_3_2:4
all.equal( colSums( vcov( estResultHom )[ 14:15, ] ), 
   - vcov( estResultHom )[ 12, ] ) # b_4_2:4
# different order of weights
estResultHom2 <- quadFuncEst( yName = "qOutput", 
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms,
   homWeights = c( qVarInput = 0.2, land = 0.1, qLabor = 0.7 ) )
all.equal( coef( estResultHom ), coef( estResultHom2 ) )
all.equal( vcov( estResultHom ), vcov( estResultHom2 ) )
# different order of xNames
estResultHom3 <- quadFuncEst( yName = "qOutput", 
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   data = germanFarms, 
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2 ) )
coefOrderHom3 <- c( 1, 5, 2:4, 15, 9, 12, 14, 6:8, 10, 11, 13 )
all.equal( coef( estResultHom ), coef( estResultHom3 )[ coefOrderHom3 ],
   check.attributes = FALSE )
all.equal( vcov( estResultHom ), 
   vcov( estResultHom3 )[ coefOrderHom3, coefOrderHom3 ],
   check.attributes = FALSE )
# homogenous in all variables
estResultHom4 <- quadFuncEst( yName = "qOutput", 
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms, 
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2, time = 0 ) )
coef( estResultHom4 )
all.equal( sum( coef( estResultHom4 )[2:4] ), - coef( estResultHom4 )[5], 
   check.attributes = FALSE ) # a_1:4
all.equal( sum( coef( estResultHom4 )[6:8] ), - coef( estResultHom4 )[9], 
   check.attributes = FALSE ) # b_1_1:4
all.equal( sum( coef( estResultHom4 )[10:12] ), - coef( estResultHom4 )[7], 
   check.attributes = FALSE ) # b_2_1:4
all.equal( sum( coef( estResultHom4 )[c(11,13,14)] ), - coef( estResultHom4 )[8], 
   check.attributes = FALSE ) # b_3_1:4
all.equal( sum( coef( estResultHom4 )[c(12,14,15)] ), - coef( estResultHom4 )[9], 
   check.attributes = FALSE ) # b_4_1:4
vcov( estResultHom4 )
all.equal( rowSums( vcov( estResultHom4 )[ , 2:4 ] ), 
   - vcov( estResultHom4 )[ , 5 ] ) # a_1:4
all.equal( rowSums( vcov( estResultHom4 )[ , 6:8 ] ), 
   - vcov( estResultHom4 )[ , 9 ] ) # b_1_1:4
all.equal( rowSums( vcov( estResultHom4 )[ , 10:12 ] ), 
   - vcov( estResultHom4 )[ , 7 ] ) # b_2_1:4
all.equal( rowSums( vcov( estResultHom4 )[ , c(11,13,14) ] ), 
   - vcov( estResultHom4 )[ , 8 ] ) # b_3_1:4
all.equal( rowSums( vcov( estResultHom4 )[ , c(12,14,15) ] ), 
   - vcov( estResultHom4 )[ , 9 ] ) # b_4_1:4
all.equal( colSums( vcov( estResultHom4 )[ 2:4, ] ), 
   - vcov( estResultHom4 )[ 5, ] ) # a_1:4
all.equal( colSums( vcov( estResultHom4 )[ 6:8, ] ), 
   - vcov( estResultHom4 )[ 9, ] ) # b_1_1:4
all.equal( colSums( vcov( estResultHom4 )[ 10:12, ] ), 
   - vcov( estResultHom4 )[ 7, ] ) # b_2_1:4
all.equal( colSums( vcov( estResultHom4 )[ c(11,13,14), ] ), 
   - vcov( estResultHom4 )[ 8, ] ) # b_3_1:4
all.equal( colSums( vcov( estResultHom4 )[ c(12,14,15), ] ), 
   - vcov( estResultHom4 )[ 9, ] ) # b_4_1:4

## computing fitted values of quadratic functions with homogeneity imposed
fitted <- quadFuncCalc(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom ), data = germanFarms, 
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2 ) )
all.equal( estResultHom$fitted, fitted )
fitted <- quadFuncCalc(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom ), data = germanFarms, 
   homWeights = c( qVarInput = 0.2, land = 0.1, qLabor = 0.7 ) )
all.equal( estResultHom$fitted, fitted )
fitted <- quadFuncCalc(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom2 ), data = germanFarms, 
   homWeights = c( qVarInput = 0.2, land = 0.1, qLabor = 0.7 ) )
all.equal( estResultHom2$fitted, fitted )
fitted <- quadFuncCalc(
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   coef = coef( estResultHom3 ), data = germanFarms,
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2 ) )
all.equal( estResultHom3$fitted, fitted )
fitted <- quadFuncCalc(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom4 ), data = germanFarms,
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2, time = 0 ) )
all.equal( estResultHom4$fitted, fitted )

## derivatives of quadratic functions with homogeneity imposed
estResultHomDeriv <- quadFuncDeriv(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom ), data = germanFarms, 
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2 ) )
print( estResultHomDeriv )
all.equal( estResultHomDeriv$qLabor * germanFarms$qLabor +
   estResultHomDeriv$land * germanFarms$land,
   - estResultHomDeriv$qVarInput * germanFarms$qVarInput )
# different order of weights (different from order used for estimation)
estResultHom1Deriv <- quadFuncDeriv(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom ), data = germanFarms, 
   homWeights = c( qVarInput = 0.2, land = 0.1, qLabor = 0.7 ) )
all.equal( estResultHomDeriv, estResultHom1Deriv )
# different order of weights (same order as used for estimation)
estResultHom2Deriv <- quadFuncDeriv(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom2 ), data = germanFarms, 
   homWeights = c( qVarInput = 0.2, land = 0.1, qLabor = 0.7 ) )
all.equal( estResultHomDeriv, estResultHom2Deriv )
# different order of independent variables
estResultHom3Deriv <- quadFuncDeriv(
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   coef = coef( estResultHom3 ), data = germanFarms,
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2 ) )
all.equal( estResultHomDeriv, 
   estResultHom3Deriv[ , c( 4, 1, 2, 3 ) ] )
# homogenous in all independent variables
estResultHom4Deriv <- quadFuncDeriv(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom4 ), data = germanFarms,
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2, time = 0 ) )
print( estResultHom4Deriv )
all.equal( estResultHom4Deriv$qLabor * germanFarms$qLabor +
   estResultHom4Deriv$land * germanFarms$land +
   estResultHom4Deriv$qVarInput * germanFarms$qVarInput,
   - estResultHom4Deriv$time * germanFarms$time )

## elasticities of quadratic functions with homogeneity imposed
estResultHomEla <- quadFuncEla(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom ), data = germanFarms,
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2 ) )
all.equal( estResultHomEla, elas( estResultHom ) )
print( estResultHomEla )
all.equal( estResultHomEla$qLabor + estResultHomEla$land,
   - estResultHomEla$qVarInput )
estResultHomElaObs <- quadFuncEla(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom ), data = germanFarms, yName = "qOutput",
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2 ) )
all.equal( estResultHomElaObs, elas( estResultHom, yObs = TRUE ) )
print( estResultHomElaObs )
all.equal( estResultHomElaObs$qLabor + estResultHomElaObs$land,
   - estResultHomElaObs$qVarInput )
max( abs( estResultHomEla - estResultHomElaObs ) )
# only at mean values
estResultHomElaMeanFit <- elas( estResultHom, data = germanFarmsMeans )
print( estResultHomElaMeanFit )
estResultHomElaMeanObs <- elas( estResultHom, data = germanFarmsMeans,
   yObs = TRUE )
print( estResultHomElaMeanObs )
print( estResultHomElaMeanFit - estResultHomElaMeanObs )
# different order of weights (different from order used for estimation)
estResultHom1Ela <- quadFuncEla(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom ), data = germanFarms,
   homWeights = c( qVarInput = 0.2, land = 0.1, qLabor = 0.7 ) )
all.equal( estResultHomEla, estResultHom1Ela )
estResultHom1ElaObs <- quadFuncEla(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom ), data = germanFarms, yName = "qOutput",
   homWeights = c( qVarInput = 0.2, land = 0.1, qLabor = 0.7 ) )
all.equal( estResultHomElaObs, estResultHom1ElaObs )
# different order of weights (same order as used for estimation)
estResultHom2Ela <- quadFuncEla(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom2 ), data = germanFarms,
   homWeights = c( qVarInput = 0.2, land = 0.1, qLabor = 0.7 ) )
all.equal( estResultHom2Ela, elas( estResultHom2 ) )
all.equal( estResultHomEla, estResultHom2Ela )
estResultHom2ElaObs <- quadFuncEla(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom2 ), data = germanFarms, yName = "qOutput",
   homWeights = c( qVarInput = 0.2, land = 0.1, qLabor = 0.7 ) )
all.equal( estResultHom2ElaObs, elas( estResultHom2, yObs = TRUE ) )
all.equal( estResultHomElaObs, estResultHom2ElaObs )
# different order of independent variables
estResultHom3Ela <- quadFuncEla(
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   coef = coef( estResultHom3 ), data = germanFarms,
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2 ) )
all.equal( estResultHom3Ela, elas( estResultHom3 ) )
all.equal( estResultHomEla, estResultHom3Ela[ , c( 4, 1, 2, 3 ) ] )
estResultHom3ElaObs <- quadFuncEla(
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   coef = coef( estResultHom3 ), data = germanFarms, yName = "qOutput",
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2 ) )
all.equal( estResultHom3ElaObs, elas( estResultHom3, yObs = TRUE ) )
all.equal( estResultHomElaObs, estResultHom3ElaObs[ , c( 4, 1, 2, 3 ) ] )
# homogenous in all independent variables
estResultHom4Ela <- quadFuncEla(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom4 ), data = germanFarms,
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2, time = 0 ) )
all.equal( estResultHom4Ela, elas( estResultHom4 ) )
print( estResultHom4Ela )
all.equal( estResultHom4Ela$qLabor + estResultHom4Ela$land +
   estResultHom4Ela$qVarInput, - estResultHom4Ela$time )
estResultHom4ElaObs <- quadFuncEla(
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHom4 ), data = germanFarms, yName = "qOutput",
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2, time = 0 ) )
all.equal( estResultHom4ElaObs, elas( estResultHom4, yObs = TRUE ) )
print( estResultHom4ElaObs )
all.equal( estResultHom4ElaObs$qLabor + estResultHom4ElaObs$land +
   estResultHom4ElaObs$qVarInput, - estResultHom4ElaObs$time )
max( abs( estResultHom4Ela - estResultHom4ElaObs ) )


################ panel data #####################
data( "GrunfeldGreene", package = "systemfit" )
ggData <- plm.data( GrunfeldGreene, c( "firm", "year" ) )
# fixed effects
ggResult <- quadFuncEst( "invest", c( "value", "capital" ), ggData )
coef( ggResult )
print( ggResult )
# random effects
ggResultRan <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   model = "random", random.method = "amemiya" )
coef( ggResultRan )
print( ggResultRan )

## panel data with a shifter
ggData$yearInt <- as.integer( as.character( ggData$year ) )
ggData$tech <- exp( ggData$yearInt - min( ggData$yearInt ) )
# fixed effects
ggResShifter <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "tech" )
coef( ggResShifter )
printQuadFuncEst( ggResShifter )
# random effects
ggResShifterRan <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "tech", model = "random", random.method = "amemiya" )
coef( ggResShifterRan )
printQuadFuncEst( ggResShifterRan )

## panel data with a logical variable as shifter
ggData$war <- ggData$yearInt >= 1939 & ggData$yearInt <= 1945
# fixed effects
ggResShifterLogi <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "war" )
coef( ggResShifterLogi )
printQuadFuncEst( ggResShifterLogi )
# random effects
ggResShifterLogiRan <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "war", model = "random", random.method = "amemiya" )
coef( ggResShifterLogiRan )
printQuadFuncEst( ggResShifterLogiRan )

## panel data with a factor as shifter
ggData$decade <- as.factor( ifelse( ggData$yearInt <= 1939, "30s",
   ifelse( ggData$yearInt <= 1949, "40s", "50s" ) ) )
# fixed effects
ggResShifterFac <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "decade" )
coef( ggResShifterFac )
printQuadFuncEst( ggResShifterFac )
# random effects
ggResShifterFacRan <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "decade", model = "random", random.method = "amemiya" )
coef( ggResShifterFacRan )
printQuadFuncEst( ggResShifterFacRan )

## linear estimations with panel data
# fixed effects
ggResultLin <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   linear = TRUE )
coef( ggResultLin )
vcov( ggResultLin )
print( ggResultLin )
# random effects
ggResultLinRan <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   linear = TRUE, model = "random", random.method = "amemiya" )
coef( ggResultLinRan )
vcov( ggResultLinRan )
print( ggResultLinRan )

## compute partial derivatives of linear estimation results with panel data
# fixed effects
margProducts <- quadFuncDeriv( c( "value", "capital" ),
   data = ggData, coef = coef( ggResultLin ), coefCov = vcov( ggResultLin ) )
sapply( margProducts, sd )
all.equal( margProducts[1,], coef( ggResultLin )[2:3],
   check.attributes = FALSE )
all.equal( attributes( margProducts )$variance[1,], diag( vcov( ggResultLin ) )[2:3], 
   check.attributes = FALSE )
# random effects
margProducts <- quadFuncDeriv( c( "value", "capital" ),
   data = ggData, coef = coef( ggResultLinRan ), coefCov = vcov( ggResultLinRan ) )
sapply( margProducts, sd )
all.equal( margProducts[1,], coef( ggResultLinRan )[2:3],
   check.attributes = FALSE )
all.equal( attributes( margProducts )$variance[1,], diag( vcov( ggResultLinRan ) )[2:3], 
   check.attributes = FALSE )

## compute elasticities using results of estimations with panel data
ggElaObs <- quadFuncEla( xNames = c( "value", "capital" ), 
   data = ggData, coef = coef( ggResult ), yName = "invest" )
all.equal( ggElaObs, elas( ggResult, yObs = TRUE ) )
print( ggElaObs )
ggElaObsRan <- quadFuncEla( xNames = c( "value", "capital" ), 
   data = ggData, coef = coef( ggResultRan ), yName = "invest" )
all.equal( ggElaObsRan, elas( ggResultRan, yObs = TRUE ) )
print( ggElaObsRan )
# with shifter variables
ggShifterElaObs <- quadFuncEla( xNames = c( "value", "capital" ), 
   data = ggData, coef = coef( ggResShifter ), yName = "invest" )
all.equal( ggShifterElaObs, elas( ggResShifter, yObs = TRUE ) )
print( ggShifterElaObs )
ggShifterElaObs2 <- quadFuncEla( xNames = c( "value", "capital" ), 
   data = ggData, coef = coef( ggResShifter ), yName = "invest",
   shifterNames = "tech" )
all.equal( ggShifterElaObs, ggShifterElaObs2 )
ggShifterElaObsRan <- quadFuncEla( xNames = c( "value", "capital" ), 
   data = ggData, coef = coef( ggResShifterRan ), yName = "invest" )
all.equal( ggShifterElaObsRan, elas( ggResShifterRan, yObs = TRUE ) )
print( ggShifterElaObsRan )

## imposing homogeneity on linear functions with panel data
ggResultLinHom <- quadFuncEst( "invest", 
   xNames = c( "value", "capital" ), data = ggData,
   linear = TRUE, homWeights = c( value = 0.3, capital = 0.7 ) )
coef( ggResultLinHom )
all.equal( coef( ggResultLinHom )[2], - coef( ggResultLinHom )[3],
   check.attributes = FALSE )
vcov( ggResultLinHom )
all.equal( vcov( ggResultLinHom )[ -1, 2 ], 
   - vcov( ggResultLinHom )[ -1, 3 ] )
all.equal( vcov( ggResultLinHom )[ 2, -1 ], 
   - vcov( ggResultLinHom )[ 3, -1 ] )

## imposing homogeneity on quadratic functions with panel data
ggResultHom <- quadFuncEst( "invest", 
   xNames = c( "value", "capital" ), data = ggData,
   homWeights = c( value = 0.3, capital = 0.7 ) )
coef( ggResultHom )
all.equal( coef( ggResultHom )[2], - coef( ggResultHom )[3],
   check.attributes = FALSE ) # a_1:2
all.equal( coef( ggResultHom )[4], - coef( ggResultHom )[5],
   check.attributes = FALSE ) # b_1_1:2
all.equal( coef( ggResultHom )[5], - coef( ggResultHom )[6],
   check.attributes = FALSE ) # b_2_1:2
vcov( ggResultHom )
all.equal( vcov( ggResultHom )[ -1, 2 ], 
   - vcov( ggResultHom )[ -1, 3 ] ) # a_1:2
all.equal( vcov( ggResultHom )[ -1, 4 ], 
   - vcov( ggResultHom )[ -1, 5 ] ) # b_1_1:2
all.equal( vcov( ggResultHom )[ -1, 5 ], 
   - vcov( ggResultHom )[ -1, 6 ] ) # b_2_1:2
all.equal( vcov( ggResultHom )[ 2, -1 ], 
   - vcov( ggResultHom )[ 3, -1 ] ) # a_1:2
all.equal( vcov( ggResultHom )[ 4, -1 ], 
   - vcov( ggResultHom )[ 5, -1 ] ) # b_1_1:2
all.equal( vcov( ggResultHom )[ 5, -1 ], 
   - vcov( ggResultHom )[ 6, -1 ] ) # b_2_1:2

########## many coefficients ##############
set.seed( 123 )
nObs <- 200
testData <- data.frame( y = rep( 3, nObs ) )
for( i in 1:15 ) {
   xName <- paste( "x", i, sep = "_" )
   testData[[ xName ]] <- rnorm( nObs )
   testData$y <- testData$y + log( i + 1 ) * testData[[ xName ]]
}
testData$y <- testData$y + 0.1 * rnorm( nObs )
testResult <- quadFuncEst( yName = "y",
   xNames = paste( "x", 1:15, sep = "_" ),
   data = testData, linear = TRUE )
coef( testResult )
rownames( vcov( testResult ) )
colnames( vcov( testResult ) )

########## scaling regressors ###########
germanFarms$landAr <- germanFarms$land * 100
germanFarms$timeCent <- germanFarms$time / 100
germanFarmsMeans <- colMeans( germanFarms[ ,
   - which( names( germanFarms ) %in% c( "year", "decade" ) ) ] )

## standard quadratic function
# estimation
estResultAr <- quadFuncEst( "qOutput",
   xNames = c( "qLabor", "landAr", "qVarInput", "time" ), data = germanFarms )
print( coef( estResult ) / coef( estResultAr ) )
# fitted values
fitted <- quadFuncCalc( xNames = c( "qLabor", "landAr", "qVarInput", "time" ),
   data = germanFarms, coef = coef( estResultAr ) )
all.equal( fitted, estResult$fitted )
# derivatives
margProductsMeanAr <- quadFuncDeriv(
   xNames = c( "qLabor", "landAr", "qVarInput", "time" ),
   data = germanFarmsMeans, coef = coef( estResultAr ),
   coefCov = vcov( estResultAr ) )
print( margProductsMean / margProductsMeanAr )
print( attributes( margProductsMean )$variance / attributes( margProductsMeanAr )$variance )
print( attributes( margProductsMean )$stdDev / attributes( margProductsMeanAr )$stdDev )
# elasticities
estElaFitAr <- elas( estResultAr )
all.equal( estElaFit, estElaFitAr, check.attributes = FALSE )
estElaObsAr <- elas( estResultAr, yObs = TRUE )
all.equal( estElaObs, estElaObsAr, check.attributes = FALSE )

## estimation with homogeneity imposed (non-homogenous variable scaled)
# estimation
estResultHomCent <- quadFuncEst( yName = "qOutput",
   xNames = c( "timeCent", "qLabor", "land", "qVarInput" ),
   data = germanFarms,
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2 ) )
print( coef( estResultHomCent ) / coef( estResultHom ) )
# fitted values
fitted <- quadFuncCalc(
   xNames = c( "timeCent", "qLabor", "land", "qVarInput" ),
   data = germanFarms, coef( estResultHomCent ),
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2 ) )
all.equal( fitted, estResultHom$fitted )
# derivatives
estResultHomDerivCent <- quadFuncDeriv(
   xNames = c( "timeCent", "qLabor", "land", "qVarInput" ),
   coef = coef( estResultHomCent ), data = germanFarms,
   homWeights = c( qLabor = 0.7, land = 0.1, qVarInput = 0.2 ) )
print( estResultHomDerivCent / estResultHomDeriv )
# elasticities
estResultHomElaCent <- elas( estResultHomCent )
all.equal( estResultHomEla, estResultHomElaCent, check.attributes = FALSE )
estResultHomElaObsCent <- elas( estResultHomCent, yObs = TRUE )
all.equal( estResultHomElaObs, estResultHomElaObsCent, check.attributes = FALSE )
