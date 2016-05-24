library(norm2)

## run EM on fake data with no missing values
set.seed(1234)
simdata <- data.frame(
   Y1=rnorm(6), Y2=rnorm(6), Y3=rnorm(6), X1=rnorm(6) )
emResult <- emNorm( cbind(Y1,Y2,Y3) ~ X1, data=simdata )
print( summary( emResult ) )

## impose missing values and run again
simdata$Y1[3] <- simdata$Y2[4] <- simdata$Y3[4] <- NA
emResult <- emNorm( cbind(Y1,Y2,Y3) ~ X1, data=simdata )
print( summary( emResult ) )

## redundant Y-variable
simdata$Y3 <- simdata$Y1 + simdata$Y2
emResult <- emNorm( cbind(Y1,Y2,Y3) ~ X1, data=simdata )

## redundant X-variable
set.seed(987)
simdata <- data.frame(
   Y1=rnorm(10), Y2=rnorm(10), Y3=rnorm(10), X1=rnorm(10) )
simdata$X2 <- simdata$X1 + 4
emResult <- emNorm( cbind(Y1,Y2,Y3) ~ X1 + X2, data=simdata )

## run EM for cholesterol data
data(cholesterol)
emResult <- emNorm(cholesterol)
print( summary(emResult) )

## re-run using formula notation
emResult <- emNorm( cbind(Y1,Y2,Y3) ~ 1, data=cholesterol)
print( summary(emResult) )

## re-run using Y1, Y2 as predictors
emResult <- emNorm( Y3 ~ Y1 + Y2, data=cholesterol)
print( summary(emResult) )

## find starting values for trivariate model using lm
tmp <- lm( cbind(Y1,Y2,Y3) ~ 1, data=cholesterol )
startval <- list(
   beta = tmp$coef,
   sigma = t(tmp$res) %*% tmp$res / tmp$df.res
   )
emResult <- emNorm( cbind(Y1,Y2,Y3) ~ 1, data=cholesterol,
   starting.values = startval )
print( summary(emResult) )

## run EM for marijuana data with ridge prior
data(marijuana)
emResult <- emNorm(marijuana, prior="ridge", prior.df=0.5)
print( summary(emResult) )

## run EM on flas data with uniform prior
data(flas)
emResult <- emNorm( flas )

## run EM on flas data with ridge prior
emResult <- emNorm( flas, prior="ridge", prior.df=0.5 )

## continue until convergence
emResult <- emNorm( emResult, iter.max=5000 )
print( summary( emResult, show.params=FALSE ) )

## treat completely observed variables as predictors;
## uniform prior
emResult <- emNorm(
   cbind(AGE,PRI,SEX,FLAS,SATV,SATM,ENG,HGPA,CGPA,GRD) ~ 
      LAN2 + LAN3 + LAN4 + MLAT, data=flas )
print( summary( emResult, show.params=FALSE ) )

## treat completely observed variables as predictors;
## ridge prior
emResult <- emNorm(
   cbind(AGE,PRI,SEX,FLAS,SATV,SATM,ENG,HGPA,CGPA,GRD) ~ 
      LAN2 + LAN3 + LAN4 + MLAT, data=flas,
   prior="ridge", prior.df=0.5 )
print( summary( emResult, show.params=FALSE ) )

