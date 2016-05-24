library( censReg )
library( plm )

options( digits = 5 )

printAll <- function( x ) {
   for( n in names( x ) ) {
      cat( "$", n, "\n", sep = "" )
      if( n %in% c( "estimate", "hessian", "gradientObs" ) ) {
         print( round( x[[ n ]], 2 ) )
      } else if( n %in% c( "gradient" ) ) {
         print( x[[ n ]], digits = 1 )
      } else if( ! n %in% c( "last.step" ) ) {
         print( x[[ n ]] )
      }
      cat( "\n" )
   }
   cat( "class\n" )
   print( class( x ) )
}

nId <- 15
nTime <- 4

set.seed( 123 )
pData <- data.frame(
   id = rep( paste( "F", 1:nId, sep = "_" ), each = nTime ),
   time = rep( 1980 + 1:nTime, nId ) )
pData$ui <- rep( rnorm( nId ), each = nTime )
pData$x1 <- rnorm( nId * nTime )
pData$x2 <- runif( nId * nTime )
pData$ys <- -1 + pData$ui + 2 * pData$x1 + 3 * pData$x2 + rnorm( nId * nTime )
pData$y <- ifelse( pData$ys > 0, pData$ys, 0 )
nData <- pData # save data set without information on panel structure
pData <- pdata.frame( pData, c( "id", "time" ) )


## Newton-Raphson method
randEff <- censReg( y ~ x1 + x2, data = pData )
print( randEff, digits = 1 )
print( randEff, logSigma = FALSE , digits = 1 )
print( maxLik:::summary.maxLik( randEff ), digits = 1 )
print( summary( randEff ), digits = 1 )
print( summary( randEff ), logSigma = FALSE , digits = 1 )
round( coef( randEff ), 2 )
round( coef( randEff, logSigma = FALSE ), 2 )
round( vcov( randEff ), 2 )
round( vcov( randEff, logSigma = FALSE ), 2 )
round( coef( summary( randEff ) ), 2 )
round( coef( summary( randEff ), logSigma = FALSE ), 2 )
try( margEff( randEff ) )
logLik( randEff )
nobs( randEff )
extractAIC( randEff )
printAll( randEff )


## BHHH method
randEffBhhh <- censReg( y ~ x1 + x2, data = pData, method = "BHHH" )
print( randEffBhhh, digits = 1 )
print( maxLik:::summary.maxLik( randEffBhhh ), digits = 1 )
print( summary( randEffBhhh ), digits = 1 )
printAll( randEffBhhh )


## BFGS method (optim)
randEffBfgs <- censReg( y ~ x1 + x2, data = pData, method = "BFGS" )
print( randEffBfgs, digits = 1 )
print( maxLik:::summary.maxLik( randEffBfgs ), digits = 1 )
print( summary( randEffBfgs ), digits = 1 )
printAll( randEffBfgs )


## BFGS method (R)
randEffBfgsr <- censReg( y ~ x1 + x2, data = pData, method = "BFGSR" )
print( randEffBfgsr, digits = 1 )
print( maxLik:::summary.maxLik( randEffBfgsr ), digits = 1 )
print( summary( randEffBfgsr ), digits = 1 )
printAll( randEffBfgsr )


## BHHH with starting values
randEffBhhhStart <- censReg( y ~ x1 + x2, data = pData, method = "BHHH",
   start = c( -0.4, 1.7, 2.2, -0.1, -0.01 ) )
print( randEffBhhhStart, digits = 1 )
print( summary( randEffBhhhStart ), digits = 1 )
nobs( randEffBhhhStart )


## left-censoring at 5
pData$yAdd <- pData$y + 5
randEffAdd <- censReg( yAdd ~ x1 + x2, data = pData, method = "BFGSR", left = 5 )
print( randEffAdd, digits = 1 )
print( maxLik:::summary.maxLik( randEffAdd ), digits = 1 )
print( summary( randEffAdd ), digits = 1 )
round( coef( randEffAdd ), 2 )
round( coef( randEffAdd, logSigma = FALSE ), 2 )
round( vcov( randEffAdd ), 2 )
round( vcov( randEffAdd, logSigma = FALSE ), 2 )
logLik( randEffAdd )
nobs( randEffAdd )
extractAIC( randEffAdd )
printAll( randEffAdd )


## right-censoring
pData$yNeg <- - pData$y
randEffNeg <- censReg( yNeg ~ x1 + x2, data = pData, method = "BFGSR",
   left = -Inf, right = 0 )
print( randEffNeg, digits = 1 )
print( maxLik:::summary.maxLik( randEffNeg ), digits = 1 )
print( summary( randEffNeg ), digits = 1 )
round( coef( randEffNeg ), 2 )
round( coef( randEffNeg, logSigma = FALSE ), 2 )
round( vcov( randEffNeg ), 2 )
round( vcov( randEffNeg, logSigma = FALSE ), 2 )
logLik( randEffNeg )
extractAIC( randEffNeg )
printAll( randEffNeg )


## right-censoring at -5
pData$yAddNeg <- - pData$yAdd
randEffAddNeg <- censReg( yAddNeg ~ x1 + x2, data = pData, method = "BFGSR",
   left = -Inf, right = -5 )
print( randEffAddNeg, digits = 1 )
print( maxLik:::summary.maxLik( randEffAddNeg ), digits = 1 )
print( summary( randEffAddNeg ), digits = 1 )
round( coef( randEffAddNeg ), 2 )
round( coef( randEffAddNeg, logSigma = FALSE ), 2 )
round( vcov( randEffAddNeg ), 2 )
round( vcov( randEffAddNeg, logSigma = FALSE ), 2 )
logLik( randEffAddNeg )
extractAIC( randEffAddNeg )
printAll( randEffAddNeg )


## both right and left censoring
pData$yBoth <- ifelse( pData$y < 3, pData$y, 3 )
randEffBoth <- censReg( yBoth ~ x1 + x2, data = pData, method = "BFGSR",
   left = 0, right = 3 )
print( randEffBoth, digits = 1 )
print( maxLik:::summary.maxLik( randEffBoth ), digits = 1 )
print( summary( randEffBoth ), digits = 1 )
print( summary( randEffBoth ), logSigma = FALSE , digits = 1 )
round( coef( randEffBoth ), 2 )
round( coef( randEffBoth, logSigma = FALSE ), 2 )
round( vcov( randEffBoth ), 2 )
round( vcov( randEffBoth, logSigma = FALSE ), 2 )
round( coef( summary( randEffBoth ) ), 2 )
round( coef( summary( randEffBoth ), logSigma = FALSE ), 2 )
logLik( randEffBoth )
nobs( randEffBoth )
extractAIC( randEffBoth )
printAll( randEffBoth )


## re-order observations/individuals
set.seed( 234 )
perm <- sample( nId )
nData2 <- nData
nData2$id <- NA
for( i in 1:nId ) {
   nData2$id[ nData$id == paste( "F", i, sep = "_" ) ] <-
      paste( "G", perm[ i ], sep = "_" )
}
pData2 <- pdata.frame( nData2, c( "id", "time" ) )
randEffBfgsr2 <- censReg( y ~ x1 + x2, data = pData2, method = "BFGSR" )
all.equal( randEffBfgsr2[ -c(11,13) ], randEffBfgsr[ -c(11,13) ] )
all.equal( sort( randEffBfgsr2[[ 11 ]] ), sort( randEffBfgsr[[ 11 ]] ) )

# check if the order of observations/individuals influences the likelihood values
d1c1 <- censReg( y ~ x1 + x2, data = pData, method = "BFGSR", start = coef(randEffBfgsr),
   iterlim = 0 )
all.equal( d1c1[-c(5,6,9,13,17)], randEffBfgsr[-c(5,6,9,13,17)] )
d1c1$maximum -  randEffBfgsr$maximum

d2c2 <- censReg( y ~ x1 + x2, data = pData2, method = "BFGSR", start = coef(randEffBfgsr2),
   iterlim = 0 )
all.equal( d2c2[-c(5,6,9,13,17)], randEffBfgsr2[-c(5,6,9,13,17)] )
d2c2$maximum -  randEffBfgsr2$maximum

d1c2 <- censReg( y ~ x1 + x2, data = pData, method = "BFGSR", 
   start = coef(randEffBfgsr2), iterlim = 0 )
d2c2$maximum - d1c2$maximum
d2c2$gradient - d1c2$gradient

d2c1 <- censReg( y ~ x1 + x2, data = pData2, method = "BFGSR", 
   start = coef(randEffBfgsr), iterlim = 0 )
d1c1$maximum - d2c1$maximum
d1c1$gradient - d2c1$gradient

d2c2$maximum - d2c1$maximum
d1c1$maximum - d1c2$maximum

d1cS <- censReg( y ~ x1 + x2, data = pData, method = "BFGSR", 
   start = randEffBfgsr$start, iterlim = 0 )
d2cS <- censReg( y ~ x1 + x2, data = pData2, method = "BFGSR", 
   start = randEffBfgsr$start, iterlim = 0 )
d1cS$maximum - d2cS$maximum
d1cS$gradient - d2cS$gradient


## unbalanced panel data
nDataUnb <- nData[ -c( 2, 5, 6, 8 ), ]
pDataUnb <- pdata.frame( nDataUnb, c( "id", "time" ) )
randEffBfgsrUnb <- censReg( y ~ x1 + x2, data = pDataUnb, method = "BFGSR" )
print( randEffBfgsrUnb )
print( maxLik:::summary.maxLik( randEffBfgsrUnb ), digits = 1 )
print( summary( randEffBfgsrUnb ), digits = 1 )
logLik( randEffBfgsrUnb )
extractAIC( randEffBfgsrUnb )
printAll( randEffBfgsrUnb )


## NAs in data
pDataNa <- pData
obsNa <- which( ! rownames( pData ) %in% rownames( pDataUnb ) )
pDataNa$y[ obsNa[ 1:2 ] ] <- NA
pDataNa$x1[ obsNa[ 3 ] ] <- NA
pDataNa$x2[ obsNa[ c( 1, 2, 4 ) ] ] <- NA
randEffBfgsrNa <- censReg( y ~ x1 + x2, data = pDataNa, method = "BFGSR" )
all.equal( randEffBfgsrNa[ -13 ], randEffBfgsrUnb[ -13 ] )


# returning log-likelihood contributions only (no estimations)
logLikRandEff <- censReg( y ~ x1 + x2, data = pData, start = coef( randEff ),
   logLikOnly = TRUE )
print( logLikRandEff, digits = 1 )
all.equal( sum( logLikRandEff ), c( logLik( randEff ) ) )
logLikStart <- censReg( y ~ x1 + x2, data = pData, 
   start = c( -0.4, 1.7, 2.2, -0.1, -0.01 ), logLikOnly = TRUE )
print( logLikStart )



