library( micEconSNQP )
data( germanFarms, package = "micEcon" )
options( digits = 3 )

germanFarms$qOutput   <- germanFarms$vOutput / germanFarms$pOutput
germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
germanFarms$qLabor    <- -germanFarms$qLabor
germanFarms$time      <- c( 0:19 )

pNamesT <- c( "pOutput", "pVarInput", "pLabor" )
qNamesT <- c( "qOutput", "qVarInput", "qLabor" )
fNamesT <- c( "land", "time" )

estResult <- snqProfitEst( pNamesT, qNamesT, "land", data = germanFarms )
print( estResult )
estResult$est <- summary( estResult$est )
print.default( estResult )

################ without fix inputs ##############################
estResult <- snqProfitEst( pNamesT, qNamesT, NULL, data = germanFarms )
print( estResult )
estResult$est <- summary( estResult$est )
print.default( estResult )

estResultCalc <- snqProfitCalc( pNamesT, NULL, estResult$data,
   estResult$weights, estResult$scalingFactors, estResult$coef )
print( estResultCalc )
if( max( abs( estResultCalc - estResult$fitted ) ) > 1e-5 ) {
   stop( "values from snqProfitCalc are not equal to fitted values" )
}

estResultEla <- snqProfitEla( estResult$coef$beta,
   estResult$data[ 20, pNamesT ], estResult$data[ 20, qNamesT ],
   estResult$weights, estResult$scalingFactors )
print( estResultEla )

estResultHessianDeriv <- snqProfitHessianDeriv( estResult$pMean,
   estResult$weights, nFix = 2 )
print( estResultHessianDeriv )

estResultHessian <- snqProfitHessian( estResult$coef$beta,
   estResult$data[ 20, pNamesT ], estResult$weights,
   estResult$scalingFactors )
print( estResultHessian )

########### with fix inputs, form = 0 ########################
estResult <- snqProfitEst( pNamesT, qNamesT, fNamesT, data=germanFarms )
print( estResult )
estResult$est <- summary( estResult$est )
print.default( estResult )

estResultCalc <- snqProfitCalc( pNamesT, fNamesT, estResult$data,
   estResult$weights, estResult$scalingFactors, estResult$coef )
print( estResultCalc )
if( max( abs( estResultCalc - estResult$fitted ) ) > 1e-5 ) {
   stop( "values from snqProfitCalc are not equal to fitted values" )
}

estResultEla <- snqProfitEla( estResult$coef$beta,
   estResult$data[ 20, pNamesT ], estResult$data[ 20, qNamesT ],
   estResult$weights, estResult$scalingFactors )
print( estResultEla )

estResultHessianDeriv <- snqProfitHessianDeriv( estResult$pMean,
   estResult$weights, nFix = 2 )
print( estResultHessianDeriv )

estResultHessian <- snqProfitHessian( estResult$coef$beta,
   estResult$data[ 20, pNamesT ], estResult$weights,
   estResult$scalingFactors )
print( estResultHessian )

estResultShadowprices <- snqProfitShadowPrices( pNamesT, fNamesT, estResult )
print( estResultShadowprices )

####################################################
estResult <- snqProfitEst( pNamesT, qNamesT, fNamesT, data=germanFarms, form = 1 )
print( estResult )
estResult$est <- summary( estResult$est )
print.default( estResult )

estResultCalc <- snqProfitCalc( pNamesT, fNamesT, estResult$data,
   estResult$weights, estResult$scalingFactors, estResult$coef, form = 1 )
print( estResultCalc )
if( max( abs( estResultCalc - estResult$fitted ) ) > 1e-5 ) {
   stop( "values from snqProfitCalc are not equal to fitted values" )
}

estResultEla <- snqProfitEla( estResult$coef$beta,
   estResult$data[ 20, pNamesT ], estResult$data[ 20, qNamesT ],
   estResult$weights, estResult$scalingFactors )
print( estResultEla )

estResultHessianDeriv <- snqProfitHessianDeriv( estResult$pMean,
   estResult$weights, nFix = 2, form = 1 )
print( estResultHessianDeriv )

estResultHessian <- snqProfitHessian( estResult$coef$beta,
   estResult$data[ 20, pNamesT ], estResult$weights,
   estResult$scalingFactors )
print( estResultHessian )

estResultShadowprices <- snqProfitShadowPrices( pNamesT, fNamesT, estResult )
print( estResultShadowprices )

########### without scaling prices ########################
estResultNoScale <- snqProfitEst( pNamesT, qNamesT, fNamesT, data=germanFarms,
   scalingFactors = rep( 1, 3 ) )
print( estResultNoScale )
estResultNoScale2 <- snqProfitEst( pNamesT, qNamesT, fNamesT, data=germanFarms,
   base = NULL, weights = snqProfitWeights( pNamesT, qNamesT, germanFarms ) )
print( estResultNoScale2 )
all.equal( estResultNoScale[-20], estResultNoScale2[] )

########### with manually chosen scaling factors ###################
estResultNoScale10 <- snqProfitEst( pNamesT, qNamesT, fNamesT, data=germanFarms,
   scalingFactors = rep( 10, 3 ) )
print( estResultNoScale10 )
all.equal( estResultNoScale[-c(1,2,4,5,6,7,9,14,22)], 
   estResultNoScale10[-c(1,2,4,5,6,7,9,14,22)], tol = 1e-4 )

estResultNoScale123 <- snqProfitEst( pNamesT, qNamesT, fNamesT, data=germanFarms,
   scalingFactors = c( 1, 2, 3 ) )
print( estResultNoScale123 )

