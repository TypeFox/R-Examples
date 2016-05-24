### R code from vignette source 'systemfit.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: systemfit.Rnw:124-125
###################################################
options( prompt = "R> ", ctinue = "+  " )


###################################################
### code chunk number 2: systemfit.Rnw:1442-1445 (eval = FALSE)
###################################################
## sigmaInv <- solve( residCov )
## xtOmegaInv <- crossprod( xMat, kronecker( sigmaInv, Diagonal( nObs ) ) )
## coef <- solve( xtOmegaInv %*% xMat, xtOmegaInv %*%  yVec )


###################################################
### code chunk number 3: systemfit_usage
###################################################
library( "systemfit" )
data( "Kmenta" )
attach( Kmenta )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
eqSystem <- list( demand = eqDemand, supply = eqSupply )
fitols <- systemfit( eqSystem )
print( fitols )


###################################################
### code chunk number 4: systemfit_usage
###################################################
fitsur <- systemfit( eqSystem, method = "SUR" )


###################################################
### code chunk number 5: systemfit_usage
###################################################
fit3sls  <- systemfit( eqSystem, method = "3SLS",
   inst = ~ income + farmPrice + trend )
fit3sls2 <- systemfit( eqSystem, method = "3SLS",
   inst = list( ~ farmPrice + trend, ~ income + farmPrice + trend ) )


###################################################
### code chunk number 6: systemfit_usage
###################################################
fitsur <- systemfit( eqSystem, method = "SUR", data = Kmenta )


###################################################
### code chunk number 7: systemfit_usage
###################################################
restrict <- "demand_price + supply_farmPrice = 0"
fitsurRmat <- systemfit(eqSystem, method = "SUR",
   restrict.matrix = restrict)


###################################################
### code chunk number 8: systemfit_usage
###################################################
Rmat <- matrix(0, nrow = 1, ncol = 7)
Rmat[1, 2] <- 1
Rmat[1, 6] <- 1
qvec <- c(0)
fitsurRmatNum <- systemfit(eqSystem, method = "SUR",
   restrict.matrix = Rmat, restrict.rhs = qvec)


###################################################
### code chunk number 9: systemfit_usage
###################################################
modRegMat <- matrix(0, nrow = 7, ncol = 6)
modRegMat[1:5, 1:5] <- diag(5)
modRegMat[6, 2] <- -1
modRegMat[7, 6] <- 1
fitsurRegMat <- systemfit(eqSystem, method = "SUR",
   restrict.regMat = modRegMat)


###################################################
### code chunk number 10: systemfit_usage
###################################################
all.equal( coef( fitsurRmat ), coef( fitsurRmatNum ) )
all.equal( coef( fitsurRmat ), coef( fitsurRegMat ) )


###################################################
### code chunk number 11: systemfit_usage
###################################################
fitsurit <- systemfit( eqSystem, method = "SUR", maxiter = 500 )


###################################################
### code chunk number 12: systemfit_usage
###################################################
summary( fitsur )


###################################################
### code chunk number 13: systemfit_usage
###################################################
summary( fitsur, residCov = FALSE, equations = FALSE )


###################################################
### code chunk number 14: systemfit_usage
###################################################
data( "GrunfeldGreene" )
library( "plm" )
GGPanel <- plm.data( GrunfeldGreene, c( "firm", "year" ) )
greeneSur <- systemfit( invest ~ value + capital, method = "SUR",
   data = GGPanel )


###################################################
### code chunk number 15: systemfit_usage
###################################################
greeneSurPooled <- systemfit( invest ~ value + capital, method = "SUR",
   data = GGPanel, pooled = TRUE )


###################################################
### code chunk number 16: systemfit_usage
###################################################
linearHypothesis( fitsur, Rmat, qvec, test = "FT" )

linearHypothesis( fitsur, Rmat, qvec, test = "F" )

linearHypothesis( fitsur, Rmat, qvec, test = "Chisq" )

lrtest( fitsurRmat, fitsur )


###################################################
### code chunk number 17: systemfit_usage
###################################################
fit2sls  <- systemfit( eqSystem, method = "2SLS",
   inst = ~ income + farmPrice + trend, data = Kmenta )
fit3sls <- systemfit( eqSystem, method = "3SLS",
   inst = ~ income + farmPrice + trend, data = Kmenta )
hausman.systemfit( fit2sls, fit3sls )


###################################################
### code chunk number 18: systemfit.Rnw:2116-2117
###################################################
hausmantest <- hausman.systemfit( fit2sls, fit3sls )


###################################################
### code chunk number 19: systemfit_replication
###################################################
library( "systemfit" )


###################################################
### code chunk number 20: systemfit_replication
###################################################
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )


###################################################
### code chunk number 21: systemfit_replication
###################################################
fitOls <- systemfit( system, data = Kmenta )
round( coef( summary( fitOls ) ), digits = 4 )


###################################################
### code chunk number 22: systemfit_replication
###################################################
fit2sls <- systemfit( system, method = "2SLS", inst = inst, data = Kmenta )
round( coef( summary( fit2sls ) ), digits = 4 )


###################################################
### code chunk number 23: systemfit_replication
###################################################
fit3sls <- systemfit( system, method = "3SLS", inst = inst, data = Kmenta )
round( coef( summary( fit3sls ) ), digits = 4 )


###################################################
### code chunk number 24: systemfit_replication
###################################################
fitI3sls <- systemfit( system, method = "3SLS", inst = inst, data = Kmenta,
   maxit = 250 )
round( coef( summary( fitI3sls ) ), digits = 4 )


###################################################
### code chunk number 25: systemfit_replication
###################################################
data( "KleinI" )
eqConsump  <- consump ~ corpProf + corpProfLag + wages
eqInvest   <- invest ~ corpProf + corpProfLag + capitalLag
eqPrivWage <- privWage ~ gnp + gnpLag + trend
inst <- ~ govExp + taxes + govWage + trend + capitalLag + corpProfLag + gnpLag
system <- list( Consumption = eqConsump, Investment = eqInvest,
   PrivateWages = eqPrivWage )


###################################################
### code chunk number 26: systemfit_replication
###################################################
kleinOls <- systemfit( system, data = KleinI )
round( coef( summary( kleinOls ) ), digits = 3 )


###################################################
### code chunk number 27: systemfit_replication
###################################################
klein2sls <- systemfit( system, method = "2SLS", inst = inst, data = KleinI,
   methodResidCov = "noDfCor" )
round( coef( summary( klein2sls ) ), digits = 3 )


###################################################
### code chunk number 28: systemfit_replication
###################################################
klein3sls <- systemfit( system, method = "3SLS", inst = inst, data = KleinI,
   methodResidCov = "noDfCor" )
round( coef( summary( klein3sls ) ), digits = 3 )


###################################################
### code chunk number 29: systemfit_replication
###################################################
kleinI3sls <- systemfit( system, method = "3SLS", inst = inst, data = KleinI,
   methodResidCov = "noDfCor", maxit = 500 )
round( coef( summary( kleinI3sls ) ), digits = 3 )


###################################################
### code chunk number 30: systemfit_replication
###################################################
data( "GrunfeldGreene" )
library( "plm" )
GGPanel <- plm.data( GrunfeldGreene, c( "firm", "year" ) )
formulaGrunfeld <- invest ~ value + capital


###################################################
### code chunk number 31: systemfit_replication
###################################################
greeneOls <- systemfit( formulaGrunfeld,
   data = GGPanel )
round( coef( summary( greeneOls ) ), digits = 4 )
round( sapply( greeneOls$eq, function(x){return(summary(x)$ssr/20)} ), digits = 3 )


###################################################
### code chunk number 32: systemfit_replication
###################################################
greeneOlsPooled <- systemfit( formulaGrunfeld,
   data = GGPanel, pooled = TRUE )
round( coef( summary( greeneOlsPooled$eq[[1]] ) ), digits = 4 ) #$
sum( sapply( greeneOlsPooled$eq, function(x){return(summary(x)$ssr)}) )/100


###################################################
### code chunk number 33: systemfit_replication
###################################################
greeneSur <- systemfit( formulaGrunfeld, method = "SUR",
   data = GGPanel, methodResidCov = "noDfCor" )
round( coef( summary( greeneSur ) ), digits = 4 )
round( greeneSur$residCov, digits = 3 ) #$
round( summary( greeneSur )$residCor, digits = 3 ) #$


###################################################
### code chunk number 34: systemfit_replication
###################################################
greeneSurPooled <- systemfit( formulaGrunfeld, method = "SUR",
   data = GGPanel, pooled = TRUE, methodResidCov = "noDfCor",
   residCovWeighted = TRUE )
round( coef( summary( greeneSurPooled$eq[[1]] ) ), digits = 4 ) #$

round( greeneSurPooled$residCov, digits = 3 ) #$
round( cov( residuals( greeneSurPooled ) ), digits = 3 )
round( summary( greeneSurPooled )$residCor, digits = 3 ) #$


###################################################
### code chunk number 35: systemfit_replication
###################################################
GrunfeldTheil <- subset( GrunfeldGreene,
   firm %in% c( "General Electric", "Westinghouse" ) )
GTPanel <- plm.data( GrunfeldTheil, c( "firm", "year" ) )
formulaGrunfeld <- invest ~ value + capital


###################################################
### code chunk number 36: systemfit_replication
###################################################
theilOls <- systemfit( formulaGrunfeld,
   data = GTPanel )
round( coef( summary( theilOls ) ), digits = 3 )


###################################################
### code chunk number 37: systemfit_replication
###################################################
theilSur <- systemfit( formulaGrunfeld, method = "SUR",
   data = GTPanel, methodResidCov = "noDfCor" )
round( coef( summary( theilSur ) ), digits = 3 )


###################################################
### code chunk number 38: systemfit_replication
###################################################
RMatrix <- rbind( c( 0, 1, 0, 0, -1, 0 ), c( 0, 0, 1, 0, 0, -1 ) )
linearHypothesis( theilSur, RMatrix )


###################################################
### code chunk number 39: systemfit_replication
###################################################
theilSurRestr <- systemfit(formulaGrunfeld, method = "SUR",
   data = GTPanel, methodResidCov = "noDfCor", restrict.matrix = RMatrix,
   residCovRestricted = FALSE)
round(coef(summary(theilSurRestr)), digits = 3)


###################################################
### code chunk number 40: systemfit_timings
###################################################
library( "systemfit" )
set.seed( 1 )
systemfitModel <- createSystemfitModel( nEq = 8, nReg = 10, nObs = 750 )
system.time(
   fitMatrix <- systemfit( systemfitModel$formula, method = "SUR",
      data = systemfitModel$data )
)
system.time(
   fitTrad <- systemfit( systemfitModel$formula, method = "SUR",
      data = systemfitModel$data, useMatrix = FALSE )
)
all.equal( fitMatrix, fitTrad )


###################################################
### code chunk number 41: systemfit_timings
###################################################
smallModel <- createSystemfitModel( nEq = 3, nReg = 4, nObs = 50 )
system.time(
   fitSmallMatrix <- systemfit( smallModel$formula, method = "SUR",
      data = smallModel$data, maxit = 500 )
)
system.time(
   fitSmallTrad <- systemfit( smallModel$formula, method = "SUR",
      data = smallModel$data, maxit = 500, useMatrix = FALSE )
)
all.equal( fitSmallMatrix, fitSmallTrad )


###################################################
### code chunk number 42: systemfit.Rnw:2643-2644
###################################################
options( width = 75 )


###################################################
### code chunk number 43: systemfit_sem
###################################################
library( "sem" )
library( "systemfit" )
data( "KleinI" )


###################################################
### code chunk number 44: systemfit_sem
###################################################
limlRam <- matrix(c(
   "consump  <-  corpProf",    "consump_corpProf",    NA,
   "consump  <-  corpProfLag", "consump_corpProfLag", NA,
   "consump  <-  wages",       "consump_wages",       NA,
   "invest   <-  corpProf",    "invest_corpProf",     NA,
   "invest   <-  corpProfLag", "invest_corpProfLag",  NA,
   "invest   <-  capitalLag",  "invest_capitalLag",   NA,
   "privWage <-  gnp",         "privWage_gnp",        NA,
   "privWage <-  gnpLag",      "privWage_gnpLag",     NA,
   "privWage <-  trend",       "privWage_trend",      NA,
   "consump  <-> consump",     "s11", NA,
   "privWage <-> privWage",    "s22", NA,
   "invest   <-> invest",      "s33", NA),
   ncol = 3, byrow = TRUE)
class(limlRam) <- "mod"
exogVar <- c("corpProf", "corpProfLag", "wages", "capitalLag", "trend",
   "gnp", "gnpLag")
endogVar <- c("consump", "invest", "privWage")
allVar <- c(exogVar, endogVar)

limlResult <- sem(model = limlRam, S = cov(KleinI[ -1, allVar ]),
   N = (nrow(KleinI) - 1), fixed.x = exogVar)
print(limlResult)


###################################################
### code chunk number 45: systemfit_sem
###################################################
eqConsump <- consump ~ corpProf + corpProfLag + wages
eqInvest <- invest ~ corpProf + corpProfLag + capitalLag
eqPrivWage <- privWage ~ gnp + gnpLag + trend
system <- list(consump = eqConsump, invest = eqInvest,
   privWage = eqPrivWage)
olsResult <- systemfit(system, data = KleinI)
print(olsResult)


###################################################
### code chunk number 46: systemfit_sem
###################################################
fimlRam <- rbind(limlRam,
   c("consump  <-> invest",   "s12", NA),
   c("consump  <-> privWage", "s13", NA),
   c("privWage <-> invest",   "s23", NA))
class(fimlRam) <- "mod"

fimlResult <- sem(model = fimlRam, S = cov(KleinI[ -1, allVar ]),
   N = (nrow(KleinI) - 1), fixed.x = exogVar)
print(fimlResult)


###################################################
### code chunk number 47: systemfit_sem
###################################################
surResult <- systemfit( system, method = "SUR", data = KleinI,
   methodResidCov = "noDfCor", maxit = 500 )
print( surResult )


