library( systemfit )
library( plm )
options( digits = 3 )
useMatrix <- FALSE

## Repeating the OLS and SUR estimations in Theil (1971, pp. 295, 300)
data( "GrunfeldGreene" )
GrunfeldTheil <- subset( GrunfeldGreene,
   firm %in% c( "General Electric", "Westinghouse" ) )
GrunfeldTheil <- plm.data( GrunfeldTheil, c( "firm", "year" ) )
formulaGrunfeld <- invest ~ value + capital

# OLS
theilOls <- systemfit( formulaGrunfeld, "OLS",
   data = GrunfeldTheil, useMatrix = useMatrix )
print( theilOls )
summary( theilOls )
summary( theilOls, useDfSys = TRUE, residCov = FALSE,
   equations = FALSE )
summary( theilOls, equations = FALSE )
coef( theilOls )
coef( summary(theilOls ) )
vcov( theilOls )
residuals( theilOls )
confint( theilOls )
fitted(theilOls  )
logLik( theilOls )
logLik( theilOls, residCovDiag = TRUE )
nobs( theilOls )
model.frame( theilOls )
model.matrix( theilOls )
formula( theilOls )
formula( theilOls$eq[[ 1 ]] )
terms( theilOls )
terms( theilOls$eq[[ 1 ]] )

# SUR
theilSur <- systemfit( formulaGrunfeld, "SUR",
   data = GrunfeldTheil, methodResidCov = "noDfCor", useMatrix = useMatrix )
print( theilSur )
summary( theilSur )
summary( theilSur, useDfSys = TRUE, equations = FALSE )
summary( theilSur, residCov = FALSE, equations = FALSE )
coef( theilSur )
coef( summary( theilSur ) )
vcov( theilSur )
residuals( theilSur )
confint( theilSur )
fitted( theilSur )
logLik( theilSur )
logLik( theilSur, residCovDiag = TRUE )
nobs( theilSur )
model.frame( theilSur )
model.matrix( theilSur )
formula( theilSur )
formula( theilSur$eq[[ 2 ]] )
terms( theilSur )
terms( theilSur$eq[[ 2 ]] )


## Repeating the OLS and SUR estimations in Greene (2003, pp. 351)
GrunfeldGreene <- plm.data( GrunfeldGreene, c( "firm", "year" ) )
formulaGrunfeld <- invest ~ value + capital

# OLS
greeneOls <- systemfit( formulaGrunfeld, "OLS",
   data = GrunfeldGreene, useMatrix = useMatrix )
print( greeneOls )
summary( greeneOls )
summary( greeneOls, useDfSys = TRUE, equations = FALSE )
summary( greeneOls, residCov = FALSE )
sapply( greeneOls$eq, function(x){return(summary(x)$ssr/20)} ) # sigma^2
coef( greeneOls )
coef( summary( greeneOls ) )
vcov( greeneOls )
residuals( greeneOls )
confint(greeneOls  )
fitted( greeneOls )
logLik( greeneOls )
logLik( greeneOls, residCovDiag = TRUE )
nobs( greeneOls )
model.frame( greeneOls )
model.matrix( greeneOls )
formula( greeneOls )
formula( greeneOls$eq[[ 2 ]] )
terms( greeneOls )
terms( greeneOls$eq[[ 2 ]] )

# OLS Pooled
greeneOlsPooled <- systemfit( formulaGrunfeld, "OLS",
   data = GrunfeldGreene, pooled = TRUE, useMatrix = useMatrix )
print( greeneOlsPooled )
summary( greeneOlsPooled )
summary( greeneOlsPooled, useDfSys = FALSE, residCov = FALSE )
summary( greeneOlsPooled, residCov = FALSE, equations = FALSE )
sum( sapply( greeneOlsPooled$eq, function(x){return(summary(x)$ssr)}) )/97 # sigma^2
coef( greeneOlsPooled )
coef( greeneOlsPooled, modified.regMat = TRUE )
coef( summary( greeneOlsPooled ) )
coef( summary( greeneOlsPooled ), modified.regMat = TRUE )
vcov( greeneOlsPooled )
vcov( greeneOlsPooled, modified.regMat = TRUE )
residuals( greeneOlsPooled )
confint( greeneOlsPooled )
fitted( greeneOlsPooled )
logLik( greeneOlsPooled )
logLik( greeneOlsPooled, residCovDiag = TRUE )
nobs( greeneOlsPooled )
model.frame( greeneOlsPooled )
model.matrix( greeneOlsPooled )
formula( greeneOlsPooled )
formula( greeneOlsPooled$eq[[ 1 ]] )
terms( greeneOlsPooled )
terms( greeneOlsPooled$eq[[ 1 ]] )

# SUR
greeneSur <- systemfit( formulaGrunfeld, "SUR",
   data = GrunfeldGreene, methodResidCov = "noDfCor", useMatrix = useMatrix )
print( greeneSur )
summary( greeneSur )
summary( greeneSur, useDfSys = TRUE, residCov = FALSE )
summary( greeneSur, equations = FALSE )
coef( greeneSur )
coef( summary( greeneSur ) )
vcov( greeneSur )
residuals( greeneSur )
confint( greeneSur )
fitted( greeneSur )
logLik( greeneSur )
logLik( greeneSur, residCovDiag = TRUE )
nobs( greeneSur )
model.frame( greeneSur )
model.matrix( greeneSur )
formula( greeneSur )
formula( greeneSur$eq[[ 1 ]] )
terms( greeneSur )
terms( greeneSur$eq[[ 1 ]] )

# SUR Pooled
greeneSurPooled <- systemfit( formulaGrunfeld, "SUR",
   data = GrunfeldGreene, pooled = TRUE, methodResidCov = "noDfCor",
   residCovWeighted = TRUE, useMatrix = useMatrix )
print( greeneSurPooled )
summary( greeneSurPooled )
summary( greeneSurPooled, useDfSys = FALSE, equations = FALSE )
summary( greeneSurPooled, residCov = FALSE, equations = FALSE )
coef( greeneSurPooled )
coef( greeneSurPooled, modified.regMat = TRUE )
coef( summary( greeneSurPooled ) )
coef( summary( greeneSurPooled ), modified.regMat = TRUE )
vcov( greeneSurPooled )
vcov( greeneSurPooled, modified.regMat = TRUE )
residuals( greeneSurPooled )
confint( greeneSurPooled )
fitted( greeneSurPooled )
logLik( greeneSurPooled )
logLik( greeneSurPooled, residCovDiag = TRUE )
nobs( greeneSurPooled )
model.frame( greeneSurPooled )
model.matrix( greeneSurPooled )
formula( greeneSurPooled )
formula( greeneSurPooled$eq[[ 1 ]] )
terms( greeneSurPooled )
terms( greeneSurPooled$eq[[ 1 ]] )


## **************** estfun ************************
library( "sandwich" )

estfun( theilOls )
round( colSums( estfun( theilOls ) ), digits = 7 )

estfun( theilSur )
round( colSums( estfun( theilSur ) ), digits = 7 )

estfun( greeneOls )
round( colSums( estfun( greeneOls ) ), digits = 7 )

try( estfun( greeneOlsPooled ) )

estfun( greeneSur )
round( colSums( estfun( greeneSur ) ), digits = 7 )

try( estfun( greeneSurPooled ) )


## **************** bread ************************
bread( theilOls )

bread( theilSur )

bread( greeneOls )

try( bread( greeneOlsPooled ) )

bread( greeneSur )

try( bread( greeneSurPooled ) )
