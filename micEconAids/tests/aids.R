library( "micEconAids" )
data( Blanciforti86 )
options( digits = 3 )

set <- !is.na( Blanciforti86$pFood1 )
setWo1 <- set & rownames( Blanciforti86 ) != 1947
pNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
wNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )
wMeans <- colMeans( Blanciforti86[ set, wNames ] )
pMeans <- colMeans( Blanciforti86[ set, pNames ] )
xtMean <- mean( Blanciforti86[ set, "xFood" ] )


## estimations with different price indices
# AIDS: translog
estResultTl <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], method = "IL" )
print( estResultTl )
print( summary( estResultTl ) )
nobs( estResultTl )

# LA-AIDS: Stone
estResultLaS <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "S" )
print( estResultLaS )
print( summary( estResultLaS ) )
nobs( estResultLaS )

# LA-AIDS: Stone with lagged shares
estResultLaSl <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "SL" )
print( estResultLaSl )
print( summary( estResultLaSl ) )
nobs( estResultLaSl )

# LA-AIDS: Paasche
estResultLaP <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "P" )
print( estResultLaP )
print( summary( estResultLaP ) )
nobs( estResultLaP )

# LA-AIDS: Laspeyres
estResultLaL <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "L" )
print( estResultLaL )
print( summary( estResultLaL ) )
nobs( estResultLaL )

# LA-AIDS: Laspeyres, simplified
estResultLaLs <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "Ls" )
print( estResultLaLs )
print( summary( estResultLaLs ) )
nobs( estResultLaLs )
all.equal( estResultLaL$coef$beta, estResultLaLs$coef$beta )
all.equal( estResultLaL$coef$gamma, estResultLaLs$coef$gamma )

# LA-AIDS: Tornqvist
estResultLaT <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "T" )
print( estResultLaT )
print( summary( estResultLaT ) )
nobs( estResultLaT )


cat( paste( "\nRepeating the demand analysis of Blanciforti, Green",
   "& King (1986)\n" ) )
estResultLA <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "SL" )
print( estResultLA )
nobs( estResultLA )
print( summary( estResultLA ) )
print( elas( estResultLA, method = "Ch", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultLA, method = "Ch", quantNames = wNames ) )
# imposing restrictions via restrict.regMat
estResultLATX <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], priceIndex = "SL",
   restrict.regMat = TRUE )
print( estResultLATX )
nobs( estResultLATX )
print( summary( estResultLATX ) )
print( elas( estResultLATX, method = "Ch", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultLATX, method = "Ch", quantNames = wNames ) )

## only homogeneity (no symmetry imposed)
estResultLAhom <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ set, ], priceIndex = "SL" )
print( estResultLAhom )
nobs( estResultLAhom )
print( summary( estResultLAhom ) )
print( elas( estResultLAhom, method = "Ch", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultLAhom, method = "Ch", quantNames = wNames ) )
# imposing restrictions via restrict.regMat
estResultLAhomTX <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ set, ], priceIndex = "SL",
   restrict.regMat = TRUE )
print( estResultLAhomTX )
nobs( estResultLAhomTX )
print( summary( estResultLAhomTX ) )
print( elas( estResultLAhomTX, method = "Ch", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultLAhomTX, method = "Ch", quantNames = wNames ) )

## unrestricted (no homogeneity and no symmetry imposed)
estResultLAunr <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ set, ], priceIndex = "SL" )
print( estResultLAunr )
nobs( estResultLAunr )
print( summary( estResultLAunr ) )
print( elas( estResultLAunr, method = "Ch", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultLAunr, method = "Ch", quantNames = wNames ) )
# imposing restrictions via restrict.regMat
estResultLAunrTX <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ set, ], priceIndex = "SL",
   restrict.regMat = TRUE )
print( estResultLAunrTX )
nobs( estResultLAunrTX )
print( summary( estResultLAunrTX ) )
print( elas( estResultLAunrTX, method = "Ch", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultLAunrTX, method = "Ch", quantNames = wNames ) )


#####################################################
## Estimation with demand shifters
Blanciforti86$trend <- c( 0:( nrow( Blanciforti86 ) - 1 ) )
estResultLAtrend <- aidsEst( pNames, wNames, "xFood",
   shifterNames = c( "trend" ), data = Blanciforti86[ set, ] )
print( estResultLAtrend )
nobs( estResultLAtrend )
summary( estResultLAtrend )

Blanciforti86$trend2 <- c( 0:( nrow( Blanciforti86 ) - 1 ) )^2
estResultLAtrend2 <- aidsEst( pNames, wNames, "xFood",
   shifterNames = c( "trend", "trend2" ), data = Blanciforti86[ set, ] )
print( estResultLAtrend2 )
nobs( estResultLAtrend2 )
summary( estResultLAtrend2 )


#####################################################
cat( paste( "\nRepeating the evaluation of different elasticity formulas",
   "of Green & Alston (1990): iterated AIDS\n" ) )
estResultAIDS <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ setWo1, ], method = "IL" )
print( estResultAIDS )
nobs( estResultAIDS )
print( summary( estResultAIDS ) )
print( elas( estResultAIDS, method = "AIDS", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultAIDS, method = "AIDS", quantNames = wNames ) )
# imposing restrictions via restrict.regMat
estResultAIDSTX <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ setWo1, ], method = "IL", restrict.regMat = TRUE )
print( estResultAIDSTX )
nobs( estResultAIDSTX )
print( summary( estResultAIDSTX ) )
print( elas( estResultAIDSTX, method = "AIDS", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultAIDSTX, method = "AIDS", quantNames = wNames ) )

## only homogeneity (no symmetry imposed)
estResultAIDShom <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL" )
print( estResultAIDShom )
nobs( estResultAIDShom )
print( summary( estResultAIDShom ) )
print( elas( estResultAIDShom, method = "AIDS", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultAIDShom, method = "AIDS", quantNames = wNames ) )
# imposing restrictions via restrict.regMat
estResultAIDShomTX <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL", restrict.regMat = TRUE )
print( estResultAIDShomTX )
nobs( estResultAIDShomTX )
print( summary( estResultAIDShomTX ) )
print( elas( estResultAIDShomTX, method = "AIDS", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultAIDShomTX, method = "AIDS", quantNames = wNames ) )

## unrestricted (no homogeneity and no symmetry imposed)
estResultAIDSunr <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL" )
print( estResultAIDSunr )
nobs( estResultAIDSunr )
print( summary( estResultAIDSunr ) )
print( elas( estResultAIDSunr, method = "AIDS", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultAIDSunr, method = "AIDS", quantNames = wNames ) )
# imposing restrictions via restrict.regMat
estResultAIDSunrTX <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL", restrict.regMat = TRUE )
print( estResultAIDSunrTX )
nobs( estResultAIDSunrTX )
print( summary( estResultAIDSunrTX ) )
print( elas( estResultAIDSunrTX, method = "AIDS", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultAIDSunrTX, method = "AIDS", quantNames = wNames ) )

## with NAs
estResultLaSNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86,
   priceIndex = "S" )
print( estResultLaSNa )
nobs( estResultLaSNa )
print( summary( estResultLaSNa ) )
print( elas( estResultLaSNa, method = "AIDS", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultLaSNa, method = "AIDS", quantNames = wNames ) )

estResultLaSlNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86,
   priceIndex = "SL" )
print( estResultLaSlNa )
nobs( estResultLaSlNa )
print( summary( estResultLaSlNa ) )
print( elas( estResultLaSlNa, method = "AIDS", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultLaSlNa, method = "AIDS", quantNames = wNames ) )

estResultLaLsNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86 )
print( estResultLaLsNa )
nobs( estResultLaLsNa )
print( summary( estResultLaLsNa ) )
print( elas( estResultLaLsNa, method = "AIDS", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultLaLsNa, method = "AIDS", quantNames = wNames ) )

estResultAIDSNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86, method = "IL" )
print( estResultAIDSNa )
nobs( estResultAIDSNa )
print( summary( estResultAIDSNa ) )
print( elas( estResultAIDSNa, method = "AIDS", quantNames = wNames,
   observedShares = TRUE ) )
print( elas( estResultAIDSNa, method = "AIDS", quantNames = wNames ) )


########## Elasticities ###############
cat( "\nAIDS: Elasticities\n" )
ela <- aidsElas( estResultAIDS$coef, shares = wMeans, prices = pMeans, method = "AIDS",
   coefCov = vcov( estResultAIDS ), df = df.residual( estResultAIDS ) )
print( ela )
print( summary( ela ) )
elaTX <- aidsElas( estResultAIDSTX$coef, shares = wMeans, prices = pMeans, method = "AIDS",
   coefCov = vcov( estResultAIDSTX ), df = df.residual( estResultAIDSTX ) )
print( elaTX )
print( summary( elaTX ) )
print( all.equal( ela, elaTX ) )

print( elas( estResultAIDS, observedShares = TRUE ) )
print( elas( estResultAIDS ) )
print( summary( elas( estResultAIDS, observedShares = TRUE ) ) )
print( summary( elas( estResultAIDS ) ) )

print( elas( estResultAIDSTX, observedShares = TRUE ) )
print( elas( estResultAIDSTX ) )
print( summary( elas( estResultAIDSTX, observedShares = TRUE ) ) )
print( summary( elas( estResultAIDSTX ) ) )


cat( "\nLA: Elasticity formula of non-linear AIDS\n" )
ela <- aidsElas( estResultLA$coef, shares = wMeans, prices = pMeans, method = "AIDS",
   coefCov = vcov( estResultLA ), df = df.residual( estResultLA ) )
print( ela )
print( summary( ela ) )
elaTX <- aidsElas( estResultLATX$coef, shares = wMeans, prices = pMeans, method = "AIDS",
   coefCov = vcov( estResultLATX ), df = df.residual( estResultLATX ) )
print( elaTX )
print( summary( elaTX ) )
print( all.equal( ela, elaTX ) )

print( elas( estResultLA, method = "AIDS", observedShares = TRUE ) )
print( elas( estResultLA, method = "AIDS" ) )
print( summary( elas( estResultLA, method = "AIDS", observedShares = TRUE ) ) )
print( summary( elas( estResultLA, method = "AIDS" ) ) )

print( elas( estResultLATX, method = "AIDS", observedShares = TRUE ) )
print( elas( estResultLATX, method = "AIDS" ) )
print( summary( elas( estResultLATX, method = "AIDS", observedShares = TRUE ) ) )
print( summary( elas( estResultLATX, method = "AIDS" ) ) )

elas( estResultLaS, method = "AIDS" )

elas( estResultLaSl, method = "AIDS" )

elas( estResultLaP, method = "AIDS" )

elas( estResultLaL, method = "AIDS" )

elas( estResultLaLs, method = "AIDS" )

elas( estResultLaT, method = "AIDS" )

cat( "\n********** Elasticities ***************" )
cat( "\nLA: Elasticity formula of Goddard or Chalfant\n" )
ela <- aidsElas( estResultLA$coef, shares = wMeans, method = "Go",
   coefCov = vcov( estResultLA ), df = df.residual( estResultLA ),
   priceIndex = "S" )
print( ela )
print( summary( ela ) )
ela <- aidsElas( estResultLA$coef, shares = wMeans, method = "Ch",
   coefCov = vcov( estResultLA ), df = df.residual( estResultLA ),
   priceIndex = "S" )
print( ela )
print( summary( ela ) )

print( elas( estResultLA, method = "Go", observedShares = TRUE ) )
print( elas( estResultLA, method = "Go" ) )
print( summary( elas( estResultLA, observedShares = TRUE ) ) )
print( summary( elas( estResultLA ) ) )

print( elas( estResultLATX, observedShares = TRUE ) )
print( elas( estResultLATX ) )
print( summary( elas( estResultLATX, observedShares = TRUE ) ) )
print( summary( elas( estResultLATX ) ) )

elas( estResultLaS, method = "Go" )

elas( estResultLaSl, method = "Go" )

elas( estResultLaP, method = "Go" )

elas( estResultLaL, method = "Go" )

elas( estResultLaLs, method = "Go" )

elas( estResultLaT, method = "Go" )

cat( "\nLA: Elasticity formula of Eales + Unnevehr\n" )
ela <- aidsElas( estResultLA$coef, shares = wMeans, method = "EU",
   priceIndex = "S" )
print( ela )

elas( estResultLaS, method = "EU" )

elas( estResultLaSl, method = "EU" )

elas( estResultLaP, method = "EU" )

elas( estResultLaL, method = "EU" )

elas( estResultLaLs, method = "EU" )

elas( estResultLaT, method = "EU" )

cat( "\nLA: Elasticity formula of Green + Alston\n" )
ela <- aidsElas( estResultLA$coef, shares = wMeans, prices = pMeans, method = "GA",
   priceIndex = "S" )
print( ela )

cat( "\nLA: Elasticity formula of Buse\n" )
ela <- aidsElas( estResultLA$coef, shares = wMeans, prices = pMeans, method = "B1",
   priceIndex = "S" )
print( ela )

elas( estResultLaS, method = "B1" )

elas( estResultLaSl, method = "B1" )

elas( estResultLaP, method = "B1" )

elas( estResultLaL, method = "B1" )

elas( estResultLaLs, method = "B1" )

elas( estResultLaT, method = "B1" )

cat( "\nLA: Elasticity formula of Buse (alternative formula)\n" )
ela <- aidsElas( estResultLA$coef, shares = wMeans, prices = pMeans, method = "B2",
   priceIndex = "S" )
print( ela )

elas( estResultLaS, method = "B2" )

elas( estResultLaSl, method = "B2" )

elas( estResultLaP, method = "B2" )

elas( estResultLaL, method = "B2" )

elas( estResultLaLs, method = "B2" )

elas( estResultLaT, method = "B2" )

aidsElas( coef( estResultTl ), prices = pMeans, shares = wMeans )
aidsElas( coef( estResultTl ), prices = pMeans, totExp = xtMean )

aidsElas( coef( estResultLaS ), shares = wMeans, method = "Ch",
   priceIndex = "S" )
aidsElas( coef( estResultLaS ), prices = pMeans, totExp = xtMean,
   method = "Ch", priceIndex = "S" )

aidsElas( coef( estResultLaS ), prices = pMeans, shares = wMeans,
   method = "B1", priceIndex = "S" )
aidsElas( coef( estResultLaS ), prices = pMeans, totExp = xtMean,
   method = "B1", priceIndex = "S" )


############# Price indices ##############
options( digits = 5 )
cat( "\n************** Price indices **************\n" )
cat( "\nStone index\n" )
pxS <- aidsPx( "S", pNames, shareNames = wNames, data = Blanciforti86 )
print( pxS )

cat( "\nStone index with lagged shares\n" )
pxSL <- aidsPx( "SL", pNames, shareNames = wNames, data = Blanciforti86 )
print( pxSL )

cat( "\nPaasche index\n" )
pxP <- aidsPx( "P", pNames, shareNames = wNames, data = Blanciforti86 )
print( pxP )

pxP2 <- aidsPx( "P", pNames, data = Blanciforti86, shareNames = wNames,
   base = row.names(Blanciforti86) == "1970" )
print( pxP2 )

cat( "\nLaspeyres index\n" )
pxL <- aidsPx( "L", pNames, shareNames = wNames, data = Blanciforti86 )
print( pxL )

cat( "\nLaspeyres index, simplified\n" )
pxLs <- aidsPx( "Ls", pNames, shareNames = wNames, data = Blanciforti86 )
print( pxLs )
print( pxL - pxLs )

pxLs2 <- aidsPx( "Ls", pNames, data = Blanciforti86, shareNames = wNames,
   base = c( 1:32 ) )
print( pxLs2 )

cat( "\nTornqvist index\n" )
pxT <- aidsPx( "T", pNames, shareNames = wNames, data = Blanciforti86 )
print( pxT )

pxT2 <- aidsPx( "T", pNames, data = Blanciforti86, shareNames = wNames,
   base = list( prices = rep( 100, 4 ), shares = rep( 0.25, 4 ) ) )
print( pxT2 )

cat( "\nTranslog index\n" )
pxTL <- aidsPx( "TL", pNames, shareNames = wNames, data = Blanciforti86,
   coef = c( list( alpha0 = 0 ), estResultLA$coef ) )
print( pxTL )

# Translog index with 1 demand shifter
pxTLtrend <- aidsPx( "TL", pNames, data = Blanciforti86,
   coef = c( list( alpha0 = 0 ), estResultLAtrend$coef ),
   shifterNames = c( "trend" ) )
print( pxTLtrend )

# Translog index with 2 demand shifters
pxTLtrend2 <- aidsPx( "TL", pNames, data = Blanciforti86,
   coef = c( list( alpha0 = 0 ), estResultLAtrend2$coef ),
   shifterNames = c( "trend", "trend2" ) )
print( pxTLtrend2 )


########### fitted values #################
fitted( estResultAIDS )
fitted( estResultLA )
fitted( estResultLaSNa )
fitted( estResultLaLsNa )


########### aidsCalc #################
options( digits = 3 )
fittedAIDS <- aidsCalc( pNames, "xFood", data = Blanciforti86[ -1, ],
   coef = estResultAIDS$coef )
print( fittedAIDS )
if( max( abs( fittedAIDS$shares[ !is.na( fittedAIDS$shares ) ] -
   estResultAIDS$wFitted ) ) > 1e-5 ) {
   stop( "fitted shares of AIDS are wrong" )
}
if( max( abs( fittedAIDS$quant[ !is.na( fittedAIDS$quant ) ] -
   estResultAIDS$qFitted ) ) > 1e-5 ) {
   stop( "fitted quantities of AIDS are wrong" )
}
fittedAIDSTX <- aidsCalc( pNames, "xFood", data = Blanciforti86[ -1, ],
   coef = estResultAIDSTX$coef )
print( fittedAIDSTX )
print( all.equal( fittedAIDS, fittedAIDSTX ) )

fittedTl <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = coef( estResultTl ) )
print( fittedTl )
predictedTl <- predict( estResultTl, observedShares = TRUE )
all.equal( fittedTl, predictedTl )
predictedTl2 <- predict( estResultTl )
all.equal( fittedTl, predictedTl2 )

# LA-AIDS with Stone price price index with lagged shares
fittedLA <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLA$coef, priceIndex = estResultLA$lnp )
print( fittedLA )
if( max( abs( fittedLA$shares[ -1, ] - estResultLA$wFitted[ setWo1, ] ),
      na.rm = TRUE ) > 1e-5 ) {
   stop( "fitted shares of LA-AIDS are wrong" )
}
if( max( abs( fittedLA$quant[ -1, ] - estResultLA$qFitted[ setWo1, ] ),
      na.rm = TRUE ) > 1e-5 ) {
   stop( "fitted quantities of LA-AIDS are wrong" )
}
fittedLATX <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLATX$coef, priceIndex = estResultLATX$lnp )
print( fittedLATX )
print( all.equal( fittedLA, fittedLATX ) )

# LA-AIDS with Stone price index
# obsereved shares in the Stone price index
fittedLaSNa <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaSNa$coef, priceIndex = estResultLaSNa$lnp )
print( fittedLaSNa )
all.equal( fittedLaSNa$shares, estResultLaSNa$wFitted )
all.equal( fittedLaSNa$quant, estResultLaSNa$qFitted )
predictedLaSNa <- predict( estResultLaSNa, observedShares = TRUE )
all.equal( fittedLaSNa$shares, predictedLaSNa$shares[ set, ] )
all.equal( fittedLaSNa$quant, predictedLaSNa$quant[ set, ] )
# fitted shares in the Stone price index
fittedLaSNa2 <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaSNa$coef, priceIndex = "S" )
print( fittedLaSNa2 )
B86LaSNa2 <- cbind( Blanciforti86[ set, c( pNames, "xFood" ) ],
   fittedLaSNa2$shares )
lnp <- aidsPx( "S", pNames, shareNames = c( "w1", "w2", "w3", "w4" ),
   data = B86LaSNa2 )
fittedLaSNa2b <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaSNa$coef, priceIndex = lnp )
all.equal( fittedLaSNa2, fittedLaSNa2b )
predictedLaSNa2 <- predict( estResultLaSNa )
all.equal( fittedLaSNa2$shares, predictedLaSNa2$shares[ set, ] )
all.equal( fittedLaSNa2$quant, predictedLaSNa2$quant[ set, ] )

# LA-AIDS with Stone price index with lagged shares
# observed shares in the Stone price index
fittedLaSl <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaSl$coef, priceIndex = estResultLaSl$lnp )
print( fittedLaSl )
all.equal( fittedLaSl$shares, estResultLaSl$wFitted, check.attributes = FALSE )
all.equal( fittedLaSl$quant, estResultLaSl$qFitted, check.attributes = FALSE )
predictedLaSl <- predict( estResultLaSl, observedShares = TRUE )
all.equal( fittedLaSl, predictedLaSl )
# fitted shares in the Stone price index
fittedLaSl2 <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaSl$coef, priceIndex = "SL" )
print( fittedLaSl2 )
B86LaSl2 <- cbind( Blanciforti86[ set, c( pNames, "xFood" ) ],
   fittedLaSl2$shares )
lnp <- aidsPx( "SL", pNames, shareNames = c( "w1", "w2", "w3", "w4" ),
   data = B86LaSl2 )
fittedLaSl2b <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaSl$coef, priceIndex = lnp )
all.equal( fittedLaSl2$shares[-1,], fittedLaSl2b$shares[-1,] )
all.equal( fittedLaSl2$quant[-1,], fittedLaSl2b$quant[-1,] )
predictedLaSl2 <- predict( estResultLaSl )
all.equal( fittedLaSl2, predictedLaSl2 )

# LA-AIDS with Paasche price index
# obsereved shares in the Paasche price index
fittedLaP <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaP$coef, priceIndex = estResultLaP$lnp )
print( fittedLaP )
all.equal( fittedLaP$shares, estResultLaP$wFitted, check.attributes = FALSE )
all.equal( fittedLaP$quant, estResultLaP$qFitted, check.attributes = FALSE )
predictedLaP <- predict( estResultLaP, observedShares = TRUE )
all.equal( fittedLaP, predictedLaP )
# fitted shares in the Stone price index
fittedLaP2 <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaP$coef, priceIndex = "P",
   basePrices = as.numeric( Blanciforti86[ 1, pNames ] ) )
print( fittedLaP2 )
B86LaP2 <- cbind( Blanciforti86[ set, c( pNames, "xFood" ) ],
   fittedLaP2$shares )
lnp <- aidsPx( "P", pNames, shareNames = c( "w1", "w2", "w3", "w4" ),
   data = B86LaP2 )
fittedLaP2b <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaP$coef, priceIndex = lnp )
all.equal( fittedLaP2, fittedLaP2b, check.attributes = FALSE )
predictedLaP2 <- predict( estResultLaP )
all.equal( fittedLaP2, predictedLaP2 )

# LA-AIDS with Laspeyres price index
fittedLaL <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaL$coef, priceIndex = estResultLaL$lnp )
print( fittedLaL )
all.equal( fittedLaL$shares, estResultLaL$wFitted[ set, ],
   check.attributes = FALSE )
all.equal( fittedLaL$quant, estResultLaL$qFitted[ set, ],
   check.attributes = FALSE )
predictedLaL <- predict( estResultLaL, observedShares = TRUE )
all.equal( fittedLaL, predictedLaL )
fittedLaL2 <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaL$coef, priceIndex = "L",
   basePrices = as.numeric( Blanciforti86[ 1, pNames ] ),
   baseShares = as.numeric( Blanciforti86[ 1, wNames ] ) )
all.equal( estResultLaL$wFitted, fittedLaL2$shares,
   check.attributes = FALSE )
all.equal( estResultLaL$qFitted, fittedLaL2$quant,
   check.attributes = FALSE )
predictedLaL2 <- predict( estResultLaL )
all.equal( fittedLaL2, predictedLaL2 )

# LA-AIDS with simplified Laspeyres price index and NAs
fittedLaLsNa <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaLsNa$coef, priceIndex = estResultLaLsNa$lnp )
print( fittedLaLsNa )
all.equal( fittedLaLsNa$shares, estResultLaLsNa$wFitted[ set, ],
   check.attributes = FALSE )
all.equal( fittedLaLsNa$quant, estResultLaLsNa$qFitted[ set, ],
   check.attributes = FALSE )
predictedLaLsNa <- predict( estResultLaLsNa, observedShares = TRUE )
all.equal( fittedLaLsNa$shares, predictedLaLsNa$shares[ set, ],
   check.attributes = FALSE )
all.equal( fittedLaLsNa$quant, predictedLaLsNa$quant[ set, ],
   check.attributes = FALSE )
fittedLaLsNa2 <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaLsNa$coef, priceIndex = "Ls",
   baseShares = as.numeric( Blanciforti86[ 1, wNames ] ) )
all.equal( estResultLaLsNa$wFitted, fittedLaLsNa2$shares,
   check.attributes = FALSE )
all.equal( estResultLaLsNa$qFitted, fittedLaLsNa2$quant,
   check.attributes = FALSE )
predictedLaLsNa2 <- predict( estResultLaLsNa )
all.equal( fittedLaLsNa2$shares, predictedLaLsNa2$shares[ set, ],
   check.attributes = FALSE )
all.equal( fittedLaLsNa2$quant, predictedLaLsNa2$quant[ set, ],
   check.attributes = FALSE )

# LA-AIDS with Tornqvist price index
# obsereved shares in the Paasche price index
fittedLaT <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaT$coef, priceIndex = estResultLaT$lnp )
print( fittedLaT )
all.equal( fittedLaT$shares, estResultLaT$wFitted, check.attributes = FALSE )
all.equal( fittedLaT$quant, estResultLaT$qFitted, check.attributes = FALSE )
predictedLaT <- predict( estResultLaT, observedShares = TRUE )
all.equal( fittedLaT, predictedLaT )
# fitted shares in the Tornqvist price index
fittedLaT2 <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaT$coef, priceIndex = "T",
   basePrices = as.numeric( Blanciforti86[ 1, pNames ] ),
   baseShares = as.numeric( Blanciforti86[ 1, wNames ] ) )
print( fittedLaT2 )
B86LaT2 <- cbind( Blanciforti86[ set, c( pNames, "xFood" ) ],
   fittedLaT2$shares )
lnp <- aidsPx( "T", pNames, shareNames = c( "w1", "w2", "w3", "w4" ),
   data = B86LaT2, base = list(
   prices = as.numeric( Blanciforti86[ 1, pNames ] ),
   shares = as.numeric( Blanciforti86[ 1, wNames ] ) ) )
fittedLaT2b <- aidsCalc( pNames, "xFood", data = Blanciforti86[ set, ],
   coef = estResultLaT$coef, priceIndex = lnp )
all.equal( fittedLaT2, fittedLaT2b, check.attributes = FALSE )
predictedLaT2 <- predict( estResultLaT )
all.equal( fittedLaT2, predictedLaT2 )

# prediction with different data
B86new <- Blanciforti86[ 2 * c( 1:floor( nrow( Blanciforti86 ) / 2 ) ),
   c( pNames, wNames, "xFood" ) ]
B86new$pFood1 <- B86new$pFood1 * 1.1
B86new$xFood <- B86new$xFood * 1.05

predNewAIDS <- predict( estResultTl, newdata = B86new, observedShares = TRUE )
print( predNewAIDS )
predNewAIDS2 <- predict( estResultTl, newdata = B86new )
all.equal( predNewAIDS, predNewAIDS2 )

predict( estResultLaS, newdata = B86new, observedShares = TRUE )
predict( estResultLaS, newdata = B86new )

predict( estResultLaSl, newdata = B86new, observedShares = TRUE )
predict( estResultLaSl, newdata = B86new )

predict( estResultLaP, newdata = B86new, observedShares = TRUE )
predict( estResultLaP, newdata = B86new )

predict( estResultLaL, newdata = B86new, observedShares = TRUE )
predict( estResultLaL, newdata = B86new )

predict( estResultLaLs, newdata = B86new, observedShares = TRUE )
predict( estResultLaLs, newdata = B86new )

predict( estResultLaT, newdata = B86new, observedShares = TRUE )
predict( estResultLaT, newdata = B86new )


###############  aidsUtility  #################
coefTl <- coef( estResultTl )
aidsUtility( pNames, "xFood", coef = coefTl, data = Blanciforti86[ set, ] )
coefTl$beta0 <- 2
aidsUtility( pNames, "xFood", coef = coefTl, data = Blanciforti86[ set, ] )

aidsUtility( pNames, "xFood", coef = coef( estResultAIDS ),
   data = Blanciforti86[ setWo1, ] )

aidsUtility( pNames, "xFood", coef = coef( estResultAIDShom ),
   data = Blanciforti86[ setWo1, ] )

aidsUtility( pNames, "xFood", coef = coef( estResultAIDSunr ),
   data = Blanciforti86[ setWo1, ] )

aidsUtility( pNames, "xFood", coef = coef( estResultAIDSNa ),
   data = Blanciforti86[ setWo1, ] )


###############  aidsUtilityDeriv  #################
coefTl <- coef( estResultTl )
aidsUtilityDeriv( pNames, "xFood", coef = coefTl,
   data = Blanciforti86[ set, ] )
aidsUtilityDeriv( pNames, "xFood", coef = coefTl,
   data = Blanciforti86[ set, ], rel = TRUE )

coefTl$beta0 <- 2
aidsUtilityDeriv( pNames, "xFood", coef = coefTl,
   data = Blanciforti86[ set, ] )
aidsUtilityDeriv( pNames, "xFood", coef = coefTl,
   data = Blanciforti86[ set, ], rel = TRUE )

aidsUtilityDeriv( pNames, "xFood", coef = coef( estResultAIDS ),
   data = Blanciforti86[ setWo1, ] )
aidsUtilityDeriv( pNames, "xFood", coef = coef( estResultAIDS ),
   data = Blanciforti86[ setWo1, ], rel = TRUE )

aidsUtilityDeriv( pNames, "xFood", coef = coef( estResultAIDShom ),
   data = Blanciforti86[ setWo1, ] )
aidsUtilityDeriv( pNames, "xFood", coef = coef( estResultAIDShom ),
   data = Blanciforti86[ setWo1, ], rel = TRUE )

aidsUtilityDeriv( pNames, "xFood", coef = coef( estResultAIDSunr ),
   data = Blanciforti86[ setWo1, ] )
aidsUtilityDeriv( pNames, "xFood", coef = coef( estResultAIDSunr ),
   data = Blanciforti86[ setWo1, ], rel = TRUE )

aidsUtilityDeriv( pNames, "xFood", coef = coef( estResultAIDSNa ),
   data = Blanciforti86[ setWo1, ] )
aidsUtilityDeriv( pNames, "xFood", coef = coef( estResultAIDSNa ),
   data = Blanciforti86[ setWo1, ], rel = TRUE )


####### monotonicity ###################
# AIDS
monoAids <- aidsMono( pNames, "xFood", coef = coef( estResultAIDS ),
      data = Blanciforti86[ set, ] )
print( monoAids )
class( monoAids ) <- NULL
print( monoAids )

# LA-AIDS with Tornqvist price index
# with fitted expenditure shares in the price index
monoLaT <- aidsMono( pNames, "xFood", coef = coef( estResultLaT ),
   data = Blanciforti86[ set, ], priceIndex = "T",
   basePrices = estResultLaT$basePrices,
   baseShares = estResultLaT$baseShares )
print( monoLaT )
class( monoLaT ) <- NULL
print( monoLaT )

# LA-AIDS with Tornqvist price index
# with observed expenditure shares in the price index
monoLaT2 <- aidsMono( pNames, "xFood", coef = coef( estResultLaT ),
   data = Blanciforti86[ set, ], priceIndex = estResultLaT$lnp )
print( monoLaT2 )
class( monoLaT2 ) <- NULL
print( monoLaT2 )


####### conscavity ###################
# AIDS with fitted shares
concavAids <- aidsConcav( pNames, "xFood", coef = coef( estResultAIDS ),
   data = Blanciforti86[ set, ] )
print( concavAids )
class( concavAids ) <- NULL
print( concavAids )

# AIDS with observed shares
concavAids2 <- aidsConcav( pNames, "xFood", coef = coef( estResultAIDS ),
   data = Blanciforti86[ set, ], shareNames = wNames )
print( concavAids2 )
class( concavAids2 ) <- NULL
print( concavAids2 )

# LA-AIDS
estResultLaSCoef <- coef( estResultLaS )
estResultLaSCoef$alpha0 <- 1
concavLaS <- aidsConcav( pNames, "xFood", coef = estResultLaSCoef,
   data = Blanciforti86[ set, ] )
print( concavLaS )
class( concavLaS ) <- NULL
print( concavLaS )


####### consistency ###################
# with observed expenditure shares
consist <- aidsConsist( pNames, "xFood", Blanciforti86[ set, ],
   coef = estResultAIDS$coef, shareNames = wNames )
print( consist )
class( consist ) <- NULL
class( consist$mono ) <- NULL
class( consist$concav ) <- NULL
print( consist )

# with fitted expenditure shares
consistFitted <- aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = estResultAIDS$coef )
print( consistFitted )
class( consistFitted ) <- NULL
class( consistFitted$mono ) <- NULL
class( consistFitted$concav ) <- NULL
print( consistFitted )

aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLaS ),
   priceIndex = "S" )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLaS ),
   priceIndex = estResultLaS$lnp )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLaP ),
   priceIndex = "P", basePrices = estResultLaP$basePrices )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLaP ),
   priceIndex = estResultLaP$lnp )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLaL ),
   priceIndex = "L", basePrices = estResultLaL$basePrices,
   baseShares = estResultLaL$baseShares )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLaT ),
   priceIndex = "T", basePrices = estResultLaT$basePrices,
   baseShares = estResultLaT$baseShares )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLaT ),
   priceIndex = "T", basePrices = estResultLaT$basePrices,
   baseShares = estResultLaT$baseShares,
   shareNames = wNames )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLA ),
   priceIndex = "SL" )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLATX ),
   priceIndex = "SL" )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLAhom ),
   priceIndex = "SL" )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLAhom ),
   priceIndex = estResultLAhom$lnp, shareNames = wNames )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLAhomTX ),
   priceIndex = "SL" )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLAunr ),
   priceIndex = "SL" )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLAunrTX ),
   priceIndex = "SL" )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultAIDS ) )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultAIDS ),
   shareNames = wNames )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultAIDSTX ) )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultAIDShom ) )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultAIDShomTX ) )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultAIDSunr ) )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultAIDSunrTX ) )
aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultAIDSunrTX ),
   shareNames = wNames )

## checkConsist
# AIDS
consistAids <- aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ setWo1, ], coef = coef( estResultAIDS ) )
consistAids2 <- checkConsist( estResultAIDS )
all.equal( consistAids, consistAids2 )

# AIDS with observedShares = TRUE
consistAidsSh <- aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ setWo1, ], coef = coef( estResultAIDS ),
   shareNames = wNames )
consistAidsSh2 <- checkConsist( estResultAIDS, observedShares = TRUE )
all.equal( consistAidsSh, consistAidsSh2 )

# LA-AIDS with Laspeyres index
consistLaL <- aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLaL ),
   priceIndex = "L", basePrices = estResultLaL$basePrices,
   baseShares = estResultLaL$baseShares )
consistLaL2 <- checkConsist( estResultLaL )
all.equal( consistLaL, consistLaL2 )

# LA-AIDS with Laspeyres index and observedShares = TRUE
consistLaLSh <- aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLaL ),
   priceIndex = estResultLaL$lnp, shareNames = wNames )
consistLaLSh2 <- checkConsist( estResultLaL, observedShares = TRUE )
all.equal( consistLaLSh, consistLaLSh2 )

# LA-AIDS with Stone index with lagged shares
consistLaSl <- aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLaSl ),
   priceIndex = "SL" )
consistLaSl2 <- checkConsist( estResultLaSl )
all.equal( consistLaSl, consistLaSl2 )

# LA-AIDS with Stone index with lagged shares and observedShares = TRUE
consistLaSlSh <- aidsConsist( pNames, totExpName = "xFood",
   data = Blanciforti86[ set, ], coef = coef( estResultLaSl ),
   priceIndex = estResultLaSl$lnp, shareNames = wNames )
consistLaSlSh2 <- checkConsist( estResultLaSl, observedShares = TRUE )
all.equal( consistLaSlSh, consistLaSlSh2 )


## replicating the LA-AIDS estimation of the SAS example
# loading data set
data( USMeatConsump )

# adding shifter variables for modeling seasonal effects
USMeatConsump$co1 <- cos( 1 / 2 * 3.14159 * USMeatConsump$t )
USMeatConsump$si1 <- sin( 1 / 2 * 3.14159 * USMeatConsump$t )

# Scaling prices by their means
USMeatConsump$beef_pm <- USMeatConsump$beef_p / mean( USMeatConsump$beef_p )
USMeatConsump$pork_pm <- USMeatConsump$pork_p / mean( USMeatConsump$pork_p )
USMeatConsump$chick_pm <- USMeatConsump$chick_p / mean( USMeatConsump$chick_p )
USMeatConsump$turkey_pm <- USMeatConsump$turkey_p / mean( USMeatConsump$turkey_p )

# Estimation of the model
meatModel <- aidsEst( c( "beef_pm", "pork_pm", "chick_pm", "turkey_pm" ),
   c( "beef_w", "pork_w", "chick_w", "turkey_w" ),
   "meat_exp", shifterNames = c( "co1", "si1", "t" ),
   method = "LA:S", data = USMeatConsump, maxiter=1000 )
meatModel
summary( meatModel )


## log likelihood values
logLik( estResultLA )
logLik( estResultLATX )
logLik( estResultLAhom )
logLik( estResultLAhomTX )
logLik( estResultLAunr )
logLik( estResultLAunrTX )
logLik( estResultLAtrend )
logLik( estResultLAtrend2 )
logLik( estResultAIDS )
logLik( estResultAIDSTX )
logLik( estResultAIDShom )
logLik( estResultAIDShomTX )
logLik( estResultAIDSunr )
logLik( estResultAIDSunrTX )
logLik( estResultLaSNa )
logLik( estResultLaSlNa )
logLik( estResultLaLsNa )
logLik( estResultAIDSNa )
logLik( meatModel )


## LR tests
lrtest( estResultLA, estResultLAhom, estResultLAunr, estResultLA )
lrtest( estResultLATX, estResultLAhomTX, estResultLAunrTX, estResultLATX )
lrtest( estResultLaLs, estResultLAtrend, estResultLAtrend2, estResultLaLs )
lrtest( estResultAIDSunr, estResultAIDShom, estResultAIDS, estResultAIDSunr )
lrtest( estResultAIDSunrTX, estResultAIDShomTX, estResultAIDSTX,
   estResultAIDSunrTX )


## comparing estimations results with different methods to impose restrictions
# estResultLA vs. estResultLATX
estResultLATX$call <- NULL
estResultLATX$est$call <- NULL
estResultLATX$est$restrict.regMat <- NULL
estResultLA$call <- NULL
estResultLA$est$call <- NULL
estResultLA$est$restrict.matrix <- NULL
estResultLA$est$restrict.rhs <- NULL
print( all.equal( estResultLA, estResultLATX ) )

# estResultLAhom vs. estResultLAhomTX
estResultLAhomTX$call <- NULL
estResultLAhomTX$est$call <- NULL
estResultLAhomTX$est$restrict.regMat <- NULL
estResultLAhom$call <- NULL
estResultLAhom$est$call <- NULL
estResultLAhom$est$restrict.matrix <- NULL
estResultLAhom$est$restrict.rhs <- NULL
print( all.equal( estResultLAhom, estResultLAhomTX ) )

# estResultLAunr vs. estResultLAunrTX
estResultLAunrTX$call <- NULL
estResultLAunrTX$est$call <- NULL
estResultLAunr$call <- NULL
estResultLAunr$est$call <- NULL
print( all.equal( estResultLAunr, estResultLAunrTX ) )

# estResultAIDS vs. estResultAIDSTX
estResultAIDSTX$call <- NULL
estResultAIDSTX$est$call <- NULL
estResultAIDSTX$est$restrict.regMat <- NULL
estResultAIDS$call <- NULL
estResultAIDS$est$call <- NULL
estResultAIDS$est$restrict.matrix <- NULL
estResultAIDS$est$restrict.rhs <- NULL
print( all.equal( estResultAIDS, estResultAIDSTX ) )

# estResultAIDShom vs. estResultAIDShomTX
estResultAIDShomTX$call <- NULL
estResultAIDShomTX$est$call <- NULL
estResultAIDShomTX$est$restrict.regMat <- NULL
estResultAIDShom$call <- NULL
estResultAIDShom$est$call <- NULL
estResultAIDShom$est$restrict.matrix <- NULL
estResultAIDShom$est$restrict.rhs <- NULL
print( all.equal( estResultAIDShom, estResultAIDShomTX ) )

# estResultAIDSunr vs. estResultAIDSunrTX
estResultAIDSunrTX$call <- NULL
estResultAIDSunrTX$est$call <- NULL
estResultAIDSunr$call <- NULL
estResultAIDSunr$est$call <- NULL
print( all.equal( estResultAIDSunr, estResultAIDSunrTX ) )
