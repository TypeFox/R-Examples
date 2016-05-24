library(norm2)

## run EM on fake data with missing values
set.seed(1234)
simdata <- data.frame(
   Y1=rnorm(10), Y2=rnorm(10), Y3=rnorm(10), X1=rnorm(10) )
simdata$Y3[7] <- simdata$Y2[3]<- NA
emResult <- emNorm( cbind(Y1,Y2,Y3) ~ X1, data=simdata )

## predictive mean imputation under MLE
imp <- impNorm( emResult, method="predict")
print( round( imp, 5 ) )

## random imputation under MLE
set.seed(1234)
imp <- impNorm( emResult, method="random")
print( round( imp, 5 ) )

## run EM on marijuana data with ridge prior, then impute
data(marijuana)
emResult <- emNorm( marijuana,  prior="ridge", prior.df=0.5 )  
set.seed(23)
imp <- impNorm( emResult, method="random")
print( round( imp, 5 ) )

## use cholesterol data
data(cholesterol)
emResult <- emNorm(cholesterol)
set.seed(23)
imp1 <- impNorm( emResult )
print( imp1 )

## formula notation
set.seed(23)
imp2 <- impNorm( cbind(Y1,Y2,Y3)~1, data=cholesterol,
   param=emResult$param )
print( imp2 )


## use data matrices
set.seed(23)
imp3 <- impNorm( cholesterol, param=emResult$param )
print( imp3 )

print( all( imp1 == imp2 ) )
print( all( imp1 == imp3 ) )
