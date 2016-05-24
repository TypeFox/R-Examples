library(norm2)

## run EM on fake data with missing values
set.seed(1234)
simdata <- data.frame(
   Y1=rnorm(10), Y2=rnorm(10), Y3=rnorm(10), X1=rnorm(10) )
simdata$Y3[7] <- simdata$Y2[3]<- NA
emResult <- emNorm( cbind(Y1,Y2,Y3) ~ X1, data=simdata )
print( emResult$loglik )
loglik1 <- loglikNorm( emResult )
emResult2 <- emNorm( emResult )
loglik2 <- emResult2$loglik[1]
print( loglik1 == loglik2 )


## use cholesterol data
data(cholesterol)
emResult <- emNorm(cholesterol)
print( loglikNorm( emResult ) )
print( loglikNorm( cbind(Y1,Y2,Y3)~1, data=cholesterol, 
   param=emResult$param ))
print( loglikNorm( cholesterol, param=emResult$param ))
print( logpostNorm( emResult ) )
print( logpostNorm( emResult, prior="uniform" ) )
print( logpostNorm( emResult, prior="ridge", prior.df=3 ) )
print( logpostNorm( emResult, prior="invwish", 
   prior.sscp=diag(c(2,1,2)), prior.df=3 ) )

