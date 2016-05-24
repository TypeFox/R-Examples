library( systemfit )
options( digits = 3 )

data( "Kmenta" )
useMatrix <- FALSE

demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
system <- list( demand = demand, supply = supply )
restrm <- matrix(0,1,7)  # restriction matrix "R"
restrm[1,3] <-  1
restrm[1,7] <- -1
restrict <- "demand_income - supply_trend = 0"
restr2m <- matrix(0,2,7)  # restriction matrix "R" 2
restr2m[1,3] <-  1
restr2m[1,7] <- -1
restr2m[2,2] <- -1
restr2m[2,5] <-  1
restr2q <- c( 0, 0.5 )  # restriction vector "q" 2
restrict2 <- c( "demand_income - supply_trend = 0",
   "- demand_price + supply_price = 0.5" )
tc <- matrix(0,7,6)
tc[1,1] <- 1
tc[2,2] <- 1
tc[3,3] <- 1
tc[4,4] <- 1
tc[5,5] <- 1
tc[6,6] <- 1
tc[7,3] <- 1
restr3m <- matrix(0,1,6)  # restriction matrix "R" 2
restr3m[1,2] <- -1
restr3m[1,5] <-  1
restr3q <- c( 0.5 )  # restriction vector "q" 2
restrict3 <- "- C2 + C5 = 0.5"


## *******  single-equation OLS estimations  *********************
lmDemand <- lm( demand, data = Kmenta )
lmSupply <- lm( supply, data = Kmenta )

## *************** WLS estimation ************************
fitwls1 <- systemfit( system, "WLS", data = Kmenta, useMatrix = useMatrix )
print( summary( fitwls1 ) )
all.equal( coef( fitwls1 ), c( coef( lmDemand ), coef( lmSupply ) ),
   check.attributes = FALSE )
all.equal( coef( summary( fitwls1 ) ),
   rbind( coef( summary( lmDemand ) ), coef( summary( lmSupply ) ) ),
   check.attributes = FALSE )
all.equal( vcov( fitwls1 ),
   as.matrix( bdiag( vcov( lmDemand ), vcov( lmSupply ) ) ),
   check.attributes = FALSE )

## *************** WLS estimation (EViews-like) ************************
fitwls1e <- systemfit( system, "WLS", data = Kmenta, methodResidCov = "noDfCor",
   x = TRUE, useMatrix = useMatrix )
print( summary( fitwls1e, useDfSys = TRUE ) )
all.equal( coef( fitwls1e ), c( coef( lmDemand ), coef( lmSupply ) ),
   check.attributes = FALSE )

## ************** WLS with cross-equation restriction ***************
fitwls2 <- systemfit( system, "WLS", data = Kmenta, restrict.matrix = restrm,
   x = TRUE, useMatrix = useMatrix )
print( summary( fitwls2 ) )
# the same with symbolically specified restrictions
fitwls2Sym <- systemfit( system, "WLS", data = Kmenta,
   restrict.matrix = restrict, x = TRUE,
   useMatrix = useMatrix )
all.equal( fitwls2, fitwls2Sym )

## ************** WLS with cross-equation restriction (EViews-like) *******
fitwls2e <- systemfit( system, "WLS", data = Kmenta, restrict.matrix = restrm,
   methodResidCov = "noDfCor", useMatrix = useMatrix )
print( summary( fitwls2e ) )

## ******* WLS with cross-equation restriction via restrict.regMat **********
fitwls3 <- systemfit( system,"WLS", data = Kmenta, restrict.regMat = tc,
   x = TRUE, useMatrix = useMatrix )
print( summary( fitwls3 ) )

## ******* WLS with cross-equation restriction via restrict.regMat (EViews-like) *****
fitwls3e <- systemfit( system,"WLS", data = Kmenta, restrict.regMat = tc,
   methodResidCov = "noDfCor", useMatrix = useMatrix )
print( summary( fitwls3e ) )

## ***** WLS with 2 cross-equation restrictions ***************
fitwls4 <- systemfit( system,"WLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, useMatrix = useMatrix )
print( summary( fitwls4 ) )
# the same with symbolically specified restrictions
fitwls4Sym <- systemfit( system, "WLS", data = Kmenta,
   restrict.matrix = restrict2, useMatrix = useMatrix )
all.equal( fitwls4, fitwls4Sym )

## ***** WLS with 2 cross-equation restrictions (EViews-like) **********
fitwls4e <- systemfit( system,"WLS", data = Kmenta, methodResidCov = "noDfCor",
   restrict.matrix = restr2m, restrict.rhs = restr2q,
   x = TRUE, useMatrix = useMatrix )
print( summary( fitwls4e ) )

## *********** WLS with 2 cross-equation restrictions via R and restrict.regMat ******
fitwls5 <- systemfit( system, "WLS", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc,
   x = TRUE, useMatrix = useMatrix )
print( summary( fitwls5 ) )
# the same with symbolically specified restrictions
fitwls5Sym <- systemfit( system, "WLS", data = Kmenta,
   restrict.matrix = restrict3, restrict.regMat = tc,
   x = TRUE, useMatrix = useMatrix )
all.equal( fitwls5, fitwls5Sym )

## *********** WLS with 2 cross-equation restrictions via R and restrict.regMat (EViews-like)
fitwls5e <- systemfit( system, "WLS", data = Kmenta, methodResidCov = "noDfCor",
   restrict.matrix = restr3m, restrict.rhs = restr3q, restrict.regMat = tc,
   useMatrix = useMatrix )
print( summary( fitwls5e ) )

## *************** iterated WLS estimation *********************
fitwlsi1 <- systemfit( system, "WLS", data = Kmenta,
   maxit = 100, useMatrix = useMatrix )
print( summary( fitwlsi1, useDfSys = TRUE ) )

## *************** iterated WLS estimation (EViews-like) ************
fitwlsi1e <- systemfit( system, "WLS", data = Kmenta, methodResidCov = "noDfCor",
   maxit = 100, x = TRUE, useMatrix = useMatrix )
print( summary( fitwlsi1e, useDfSys = TRUE ) )

## ****** iterated WLS with cross-equation restriction ***************
fitwlsi2 <- systemfit( system, "WLS", data = Kmenta, restrict.matrix = restrm,
   maxit = 100, x = TRUE, useMatrix = useMatrix )
print( summary( fitwlsi2 ) )

## ****** iterated WLS with cross-equation restriction (EViews-like) ********
fitwlsi2e <- systemfit( system, "WLS", data = Kmenta, restrict.matrix = restrm,
   methodResidCov = "noDfCor", maxit = 100, useMatrix = useMatrix )
print( summary( fitwlsi2e ) )

## ******* iterated WLS with cross-equation restriction via restrict.regMat **********
fitwlsi3 <- systemfit( system, "WLS", data = Kmenta, restrict.regMat = tc,
   maxit = 100, x = TRUE, useMatrix = useMatrix )
print( summary( fitwlsi3 ) )

## ******* iterated WLS with cross-equation restriction via restrict.regMat (EViews-like) ***
fitwlsi3e <- systemfit( system, "WLS", data = Kmenta, restrict.regMat = tc,
   methodResidCov = "noDfCor", maxit = 100, useMatrix = useMatrix )
print( summary( fitwlsi3e ) )
nobs( fitwlsi3e )

## ******* iterated WLS with 2 cross-equation restrictions ***********
fitwlsi4 <- systemfit( system, "WLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, maxit = 100, useMatrix = useMatrix )
print( summary( fitwlsi4 ) )

## ******* iterated WLS with 2 cross-equation restrictions (EViews-like) *****
fitwlsi4e <- systemfit( system, "WLS", data = Kmenta, methodResidCov = "noDfCor",
   restrict.matrix = restr2m, restrict.rhs = restr2q, maxit = 100,
   x = TRUE, useMatrix = useMatrix )
print( summary( fitwlsi4e ) )

## ***** iterated WLS with 2 cross-equation restrictions via R and restrict.regMat ******
fitwlsi5 <- systemfit( system, "WLS", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, maxit = 100,
   x = TRUE, useMatrix = useMatrix )
print( summary( fitwlsi5 ) )

## *** iterated WLS with 2 cross-equation restrictions via R and restrict.regMat (EViews-like)
fitwlsi5e <- systemfit( system, "WLS", data = Kmenta, methodResidCov = "noDfCor",
   restrict.matrix = restr3m, restrict.rhs = restr3q, restrict.regMat = tc,
   maxit = 100, useMatrix = useMatrix )
print( summary( fitwlsi5e ) )


## *********** estimations with a single regressor ************
fitwlsS1 <- systemfit(
   list( consump ~ price - 1, consump ~ price + trend ), "WLS",
   data = Kmenta, useMatrix = useMatrix )
print( summary( fitwlsS1 ) )
fitwlsS2 <- systemfit(
   list( consump ~ price - 1, consump ~ trend - 1 ), "WLS",
   data = Kmenta, useMatrix = useMatrix )
print( summary( fitwlsS2 ) )
fitwlsS3 <- systemfit(
   list( consump ~ trend - 1, price ~ trend - 1 ), "WLS",
   data = Kmenta, useMatrix = useMatrix )
print( summary( fitwlsS3 ) )
fitwlsS4 <- systemfit(
   list( consump ~ trend - 1, price ~ trend - 1 ), "WLS",
   data = Kmenta, restrict.matrix = matrix( c( 1, -1 ), nrow = 1 ),
   useMatrix = useMatrix )
print( summary( fitwlsS4 ) )
fitwlsS5 <- systemfit(
   list( consump ~ 1, price ~ 1 ), "WLS",
   data = Kmenta, useMatrix = useMatrix )
print( summary( fitwlsS5) )


## **************** shorter summaries **********************
print( summary( fitwls1 ), residCov = FALSE, equations = FALSE )

print( summary( fitwls2e, useDfSys = FALSE, residCov = FALSE ),
   equations = FALSE )

print( summary( fitwls3 ), residCov = FALSE )

print( summary( fitwls4e, residCov = FALSE, equations = FALSE ) )

print( summary( fitwls5, useDfSys = FALSE ), residCov = FALSE )

print( summary( fitwlsi1e, useDfSys = TRUE, equations = FALSE ) )

print( summary( fitwlsi2, equations = FALSE, residCov = FALSE ),
   residCov = TRUE )

print( summary( fitwlsi3e ), equations = FALSE, residCov = FALSE )

print( summary( fitwlsi4, equations = FALSE ), equations = TRUE )

print( summary( fitwlsi5e, useDfSys = FALSE, residCov = FALSE ) )


## ****************** residuals **************************
print( residuals( fitwls1 ) )
print( residuals( fitwls1$eq[[ 2 ]] ) )

print( residuals( fitwls2e ) )
print( residuals( fitwls2e$eq[[ 1 ]] ) )

print( residuals( fitwls3 ) )
print( residuals( fitwls3$eq[[ 2 ]] ) )

print( residuals( fitwls4e ) )
print( residuals( fitwls4e$eq[[ 1 ]] ) )

print( residuals( fitwls5 ) )
print( residuals( fitwls5$eq[[ 2 ]] ) )

print( residuals( fitwlsi1e ) )
print( residuals( fitwlsi1e$eq[[ 1 ]] ) )

print( residuals( fitwlsi2 ) )
print( residuals( fitwlsi2$eq[[ 2 ]] ) )

print( residuals( fitwlsi3e ) )
print( residuals( fitwlsi3e$eq[[ 1 ]] ) )

print( residuals( fitwlsi4 ) )
print( residuals( fitwlsi4$eq[[ 2 ]] ) )

print( residuals( fitwlsi5e ) )
print( residuals( fitwlsi5e$eq[[ 1 ]] ) )


## *************** coefficients *********************
print( round( coef( fitwls1e ), digits = 6 ) )
print( round( coef( fitwls1e$eq[[ 1 ]] ), digits = 6 ) )

print( round( coef( fitwlsi2 ), digits = 6 ) )
print( round( coef( fitwlsi2$eq[[ 2 ]] ), digits = 6 ) )

print( round( coef( fitwls3e ), digits = 6 ) )
print( round( coef( fitwls3e, modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( fitwls3e$eq[[ 1 ]] ), digits = 6 ) )

print( round( coef( fitwls4 ), digits = 6 ) )
print( round( coef( fitwls4$eq[[ 2 ]] ), digits = 6 ) )

print( round( coef( fitwlsi5 ), digits = 6 ) )
print( round( coef( fitwlsi5, modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( fitwlsi5$eq[[ 1 ]] ), digits = 6 ) )


## *************** coefficients with stats *********************
print( round( coef( summary( fitwls1e ) ), digits = 6 ) )
print( round( coef( summary( fitwls1e$eq[[ 1 ]] ) ), digits = 6 ) )

print( round( coef( summary( fitwlsi2 ) ), digits = 6 ) )
print( round( coef( summary( fitwlsi2$eq[[ 2 ]] ) ), digits = 6 ) )

print( round( coef( summary( fitwls3e ) ), digits = 6 ) )
print( round( coef( summary( fitwls3e ), modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( summary( fitwls3e$eq[[ 1 ]] ) ), digits = 6 ) )

print( round( coef( summary( fitwls4, useDfSys = FALSE ) ), digits = 6 ) )
print( round( coef( summary( fitwls4$eq[[ 2 ]], useDfSys = FALSE ) ),
   digits = 6 ) )

print( round( coef( summary( fitwlsi5, useDfSys = FALSE ) ), digits = 6 ) )
print( round( coef( summary( fitwlsi5, useDfSys = FALSE ),
   modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( summary( fitwlsi5$eq[[ 1 ]], useDfSys = FALSE ) ),
   digits = 6 ) )


## *********** variance covariance matrix of the coefficients *******
print( round( vcov( fitwls1e ), digits = 6 ) )
print( round( vcov( fitwls1e$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitwls2 ), digits = 6 ) )
print( round( vcov( fitwls2$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitwls3e ), digits = 6 ) )
print( round( vcov( fitwls3e, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fitwls3e$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitwls4 ), digits = 6 ) )
print( round( vcov( fitwls4$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitwls5 ), digits = 6 ) )
print( round( vcov( fitwls5, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fitwls5$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitwlsi1 ), digits = 6 ) )
print( round( vcov( fitwlsi1$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitwlsi2e ), digits = 6 ) )
print( round( vcov( fitwlsi2e$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitwlsi3 ), digits = 6 ) )
print( round( vcov( fitwlsi3, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fitwlsi3$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitwlsi4e ), digits = 6 ) )
print( round( vcov( fitwlsi4e$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitwlsi5e ), digits = 6 ) )
print( round( vcov( fitwlsi5e, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fitwlsi5e$eq[[ 2 ]] ), digits = 6 ) )


## *********** confidence intervals of coefficients *************
print( confint( fitwls1 ) )
print( confint( fitwls1$eq[[ 2 ]], level = 0.9 ) )

print( confint( fitwls2e, level = 0.9 ) )
print( confint( fitwls2e$eq[[ 1 ]], level = 0.99 ) )

print( confint( fitwls3, level = 0.99 ) )
print( confint( fitwls3$eq[[ 2 ]], level = 0.5 ) )

print( confint( fitwls4e, level = 0.5 ) )
print( confint( fitwls4e$eq[[ 1 ]], level = 0.25 ) )

print( confint( fitwls5, level = 0.25 ) )
print( confint( fitwls5$eq[[ 2 ]], level = 0.975 ) )

print( confint( fitwlsi1e, level = 0.975, useDfSys = TRUE ) )
print( confint( fitwlsi1e$eq[[ 1 ]], level = 0.999, useDfSys = TRUE ) )

print( confint( fitwlsi2, level = 0.999 ) )
print( confint( fitwlsi2$eq[[ 2 ]], level = 0.1 ) )

print( confint( fitwlsi3e, level = 0.1 ) )
print( confint( fitwlsi3e$eq[[ 1 ]], level = 0.01 ) )

print( confint( fitwlsi4, level = 0.01 ) )
print( confint( fitwlsi4$eq[[ 2 ]], level = 0.33 ) )

print( confint( fitwlsi5e, level = 0.33 ) )
print( confint( fitwlsi5e$eq[[ 1 ]] ) )


## *********** fitted values *************
print( fitted( fitwls1 ) )
print( fitted( fitwls1$eq[[ 2 ]] ) )

print( fitted( fitwls2e ) )
print( fitted( fitwls2e$eq[[ 1 ]] ) )

print( fitted( fitwls3 ) )
print( fitted( fitwls3$eq[[ 2 ]] ) )

print( fitted( fitwls4e ) )
print( fitted( fitwls4e$eq[[ 1 ]] ) )

print( fitted( fitwls5 ) )
print( fitted( fitwls5$eq[[ 2 ]] ) )

print( fitted( fitwlsi1e ) )
print( fitted( fitwlsi1e$eq[[ 1 ]] ) )

print( fitted( fitwlsi2 ) )
print( fitted( fitwlsi2$eq[[ 2 ]] ) )

print( fitted( fitwlsi3e ) )
print( fitted( fitwlsi3e$eq[[ 1 ]] ) )

print( fitted( fitwlsi4 ) )
print( fitted( fitwlsi4$eq[[ 2 ]] ) )

print( fitted( fitwlsi5e ) )
print( fitted( fitwlsi5e$eq[[ 1 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$consump <- NULL
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fitwls1, se.fit = TRUE, interval = "prediction" ) )
print( predict( fitwls1$eq[[ 2 ]], se.fit = TRUE, interval = "prediction" ) )

print( predict( fitwls2e, se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )
print( predict( fitwls2e$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )

print( predict( fitwls3, se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )
print( predict( fitwls3$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )

print( predict( fitwls4e, se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )
print( predict( fitwls4e$eq[[ 1 ]], se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )

print( predict( fitwls5, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )
print( predict( fitwls5$eq[[ 2 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )

print( predict( fitwlsi1e, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99, useDfSys = TRUE ) )
print( predict( fitwlsi1e$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99, useDfSys = TRUE ) )

print( predict( fitwlsi2, se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )
print( predict( fitwlsi2$eq[[ 2 ]], se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )

print( predict( fitwlsi3e, interval = "prediction", level = 0.925 ) )
print( predict( fitwlsi3e$eq[[ 1 ]], interval = "prediction", level = 0.925 ) )

print( predict( fitwlsi4, interval = "confidence", newdata = predictData ) )
print( predict( fitwlsi4$eq[[ 2 ]], interval = "confidence",
   newdata = predictData ) )

print( predict( fitwlsi5e, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.01 ) )
print( predict( fitwlsi5e$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.01 ) )

# predict just one observation
smallData <- data.frame( price = 130, income = 150, farmPrice = 120,
   trend = 25 )

print( predict( fitwls1, newdata = smallData ) )
print( predict( fitwls1$eq[[ 1 ]], newdata = smallData ) )

print( predict( fitwls2e, se.fit = TRUE, level = 0.9,
   newdata = smallData ) )
print( predict( fitwls2e$eq[[ 1 ]], se.pred = TRUE, level = 0.99,
   newdata = smallData ) )

print( predict( fitwls3, interval = "prediction", level = 0.975,
   newdata = smallData ) )
print( predict( fitwls3$eq[[ 1 ]], interval = "confidence", level = 0.8,
   newdata = smallData ) )

print( predict( fitwls4e, se.fit = TRUE, interval = "confidence",
   level = 0.999, newdata = smallData ) )
print( predict( fitwls4e$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.75, newdata = smallData ) )

print( predict( fitwls5, se.fit = TRUE, interval = "prediction",
   newdata = smallData ) )
print( predict( fitwls5$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   newdata = smallData ) )

print( predict( fitwlsi3e, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = smallData ) )
print( predict( fitwlsi3e$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, newdata = smallData ) )


## ************ correlation of predicted values ***************
print( correlation.systemfit( fitwls1, 2, 1 ) )

print( correlation.systemfit( fitwls2e, 1, 2 ) )

print( correlation.systemfit( fitwls3, 2, 1 ) )

print( correlation.systemfit( fitwls4e, 1, 2 ) )

print( correlation.systemfit( fitwls5, 2, 1 ) )

print( correlation.systemfit( fitwlsi1e, 1, 2 ) )

print( correlation.systemfit( fitwlsi2, 2, 1 ) )

print( correlation.systemfit( fitwlsi3e, 1, 2 ) )

print( correlation.systemfit( fitwlsi4, 2, 1 ) )

print( correlation.systemfit( fitwlsi5e, 1, 2 ) )


## ************ Log-Likelihood values ***************
print( logLik( fitwls1 ) )
print( logLik( fitwls1, residCovDiag = TRUE ) )
all.equal( logLik( fitwls1, residCovDiag = TRUE ),
   logLik( lmDemand ) + logLik( lmSupply ),
   check.attributes = FALSE )

print( logLik( fitwls2e ) )
print( logLik( fitwls2e, residCovDiag = TRUE ) )

print( logLik( fitwls3 ) )
print( logLik( fitwls3, residCovDiag = TRUE ) )

print( logLik( fitwls4e ) )
print( logLik( fitwls4e, residCovDiag = TRUE ) )

print( logLik( fitwls5 ) )
print( logLik( fitwls5, residCovDiag = TRUE ) )

print( logLik( fitwlsi1e ) )
print( logLik( fitwlsi1e, residCovDiag = TRUE ) )

print( logLik( fitwlsi2 ) )
print( logLik( fitwlsi2, residCovDiag = TRUE ) )

print( logLik( fitwlsi3e ) )
print( logLik( fitwlsi3e, residCovDiag = TRUE ) )

print( logLik( fitwlsi4 ) )
print( logLik( fitwlsi4, residCovDiag = TRUE ) )

print( logLik( fitwlsi5e ) )
print( logLik( fitwlsi5e, residCovDiag = TRUE ) )


## ************** F tests ****************
# testing first restriction
print( linearHypothesis( fitwls1, restrm ) )
linearHypothesis( fitwls1, restrict )

print( linearHypothesis( fitwlsi1e, restrm ) )
linearHypothesis( fitwlsi1e, restrict )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
restrictOnly2 <- "- demand_price + supply_price = 0.5"
# first restriction not imposed 
print( linearHypothesis( fitwls1e, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitwls1e, restrictOnly2 )

print( linearHypothesis( fitwlsi1, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitwlsi1, restrictOnly2 )

# first restriction imposed
print( linearHypothesis( fitwls2, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitwls2, restrictOnly2 )

print( linearHypothesis( fitwls3, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitwls3, restrictOnly2 )

print( linearHypothesis( fitwlsi2e, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitwlsi2e, restrictOnly2 )

print( linearHypothesis( fitwlsi3e, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitwlsi3e, restrictOnly2 )

# testing both of the restrictions
print( linearHypothesis( fitwls1e, restr2m, restr2q ) )
linearHypothesis( fitwls1e, restrict2 )

print( linearHypothesis( fitwlsi1, restr2m, restr2q ) )
linearHypothesis( fitwlsi1, restrict2 )


## ************** Wald tests ****************
# testing first restriction
print( linearHypothesis( fitwls1, restrm, test = "Chisq" ) )
linearHypothesis( fitwls1, restrict, test = "Chisq" )

print( linearHypothesis( fitwlsi1e, restrm, test = "Chisq" ) )
linearHypothesis( fitwlsi1e, restrict, test = "Chisq" )

# testing second restriction
# first restriction not imposed
print( linearHypothesis( fitwls1e, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitwls1e, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fitwlsi1, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitwlsi1, restrictOnly2, test = "Chisq" )

# first restriction imposed
print( linearHypothesis( fitwls2, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitwls2, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fitwls3, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitwls3, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fitwlsi2e, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitwlsi2e, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fitwlsi3e, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitwlsi3e, restrictOnly2, test = "Chisq" )

# testing both of the restrictions
print( linearHypothesis( fitwls1e, restr2m, restr2q, test = "Chisq" ) )
linearHypothesis( fitwls1e, restrict2, test = "Chisq" )

print( linearHypothesis( fitwlsi1, restr2m, restr2q, test = "Chisq" ) )
linearHypothesis( fitwlsi1, restrict2, test = "Chisq" )


## ****************** model frame **************************
print( mf <- model.frame( fitwls1 ) )
print( mf1 <- model.frame( fitwls1$eq[[ 1 ]] ) )
print( attributes( mf1 )$terms )
print( mf2 <- model.frame( fitwls1$eq[[ 2 ]] ) )
print( attributes( mf2 )$terms )

print( all.equal( mf, model.frame( fitwls2e ) ) )
print( all.equal( mf1, model.frame( fitwls2e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitwls3 ) ) )
print( all.equal( mf2, model.frame( fitwls3$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitwls4e ) ) )
print( all.equal( mf1, model.frame( fitwls4e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitwls5 ) ) )
print( all.equal( mf2, model.frame( fitwls5$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitwlsi1e ) ) )
print( all.equal( mf1, model.frame( fitwlsi1e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitwlsi2 ) ) )
print( all.equal( mf2, model.frame( fitwlsi2$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitwlsi3e ) ) )
print( all.equal( mf1, model.frame( fitwlsi3e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitwlsi4 ) ) )
print( all.equal( mf2, model.frame( fitwlsi4$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitwlsi5e ) ) )
print( all.equal( mf1, model.frame( fitwlsi5e$eq[[ 1 ]] ) ) )


## **************** model matrix ************************
# with x (returnModelMatrix) = TRUE
print( !is.null( fitwls1e$eq[[ 1 ]]$x ) )
print( mm <- model.matrix( fitwlsi1e ) )
print( mm1 <- model.matrix( fitwlsi1e$eq[[ 1 ]] ) )
print( mm2 <- model.matrix( fitwlsi1e$eq[[ 2 ]] ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fitwlsi1 ) ) )
print( all.equal( mm1, model.matrix( fitwlsi1$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwlsi1$eq[[ 2 ]] ) ) )
print( !is.null( fitwls1$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fitwls2$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fitwls2 ) ) )
print( all.equal( mm1, model.matrix( fitwls2$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwls2$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fitwls2e ) ) )
print( all.equal( mm1, model.matrix( fitwls2e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwls2e$eq[[ 2 ]] ) ) )
print( !is.null( fitwls2e$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fitwlsi3$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fitwlsi3 ) ) )
print( all.equal( mm1, model.matrix( fitwlsi3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwlsi3$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fitwlsi3e ) ) )
print( all.equal( mm1, model.matrix( fitwlsi3e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwlsi3e$eq[[ 2 ]] ) ) )
print( !is.null( fitwlsi3e$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fitwls4e$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fitwls4e ) ) )
print( all.equal( mm1, model.matrix( fitwls4e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwls4e$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fitwls4Sym ) ) )
print( all.equal( mm1, model.matrix( fitwls4Sym$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwls4Sym$eq[[ 2 ]] ) ) )
print( !is.null( fitwls4Sym$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fitwls5$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fitwls5 ) ) )
print( all.equal( mm1, model.matrix( fitwls5$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwls5$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fitwls5e ) ) )
print( all.equal( mm1, model.matrix( fitwls5e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwls5e$eq[[ 2 ]] ) ) )
print( !is.null( fitwls5e$eq[[ 1 ]]$x ) )


## **************** formulas ************************
formula( fitwls1 )
formula( fitwls1$eq[[ 2 ]] )

formula( fitwls2e )
formula( fitwls2e$eq[[ 1 ]] )

formula( fitwls3 )
formula( fitwls3$eq[[ 2 ]] )

formula( fitwls4e )
formula( fitwls4e$eq[[ 1 ]] )

formula( fitwls5 )
formula( fitwls5$eq[[ 2 ]] )

formula( fitwlsi1e )
formula( fitwlsi1e$eq[[ 1 ]] )

formula( fitwlsi2 )
formula( fitwlsi2$eq[[ 2 ]] )

formula( fitwlsi3e )
formula( fitwlsi3e$eq[[ 1 ]] )

formula( fitwlsi4 )
formula( fitwlsi4$eq[[ 2 ]] )

formula( fitwlsi5e )
formula( fitwlsi5e$eq[[ 1 ]] )


## **************** model terms *******************
terms( fitwls1 )
terms( fitwls1$eq[[ 2 ]] )

terms( fitwls2e )
terms( fitwls2e$eq[[ 1 ]] )

terms( fitwls3 )
terms( fitwls3$eq[[ 2 ]] )

terms( fitwls4e )
terms( fitwls4e$eq[[ 1 ]] )

terms( fitwls5 )
terms( fitwls5$eq[[ 2 ]] )

terms( fitwlsi1e )
terms( fitwlsi1e$eq[[ 1 ]] )

terms( fitwlsi2 )
terms( fitwlsi2$eq[[ 2 ]] )

terms( fitwlsi3e )
terms( fitwlsi3e$eq[[ 1 ]] )

terms( fitwlsi4 )
terms( fitwlsi4$eq[[ 2 ]] )

terms( fitwlsi5e )
terms( fitwlsi5e$eq[[ 1 ]] )


## **************** estfun ************************
library( "sandwich" )

estfun( fitwls1 )
round( colSums( estfun( fitwls1 ) ), digits = 7 )

estfun( fitwlsi1e )
round( colSums( estfun( fitwlsi1e ) ), digits = 7 )


## **************** bread ************************
bread( fitwls1 )

bread( fitwlsi1e )
