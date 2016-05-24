library(norm2)

## run EM on fake data with missing values
set.seed(1234)
simdata <- data.frame(
   Y1=rnorm(10), Y2=rnorm(10), Y3=rnorm(10), X1=rnorm(10) )
simdata$Y3[7] <- simdata$Y2[3]<- NA
emResult <- emNorm( cbind(Y1,Y2,Y3) ~ X1, data=simdata )

## run mcmc
set.seed(567)
mcmcResult1 <- mcmcNorm( emResult )
acf.worst1 <- acf( mcmcResult1$series.worst, plot=FALSE )
print( round( acf.worst1$acf, 5 ) )

## re-run identical mcmc series using formula notation
set.seed(567)
mcmcResult2 <- mcmcNorm( cbind(Y1,Y2,Y3) ~ X1, data=simdata, 
   starting.values=emResult$param,
   save.worst.series=TRUE, worst.linear.coef=emResult$worst.linear.coef )
acf.worst2 <- acf( mcmcResult2$series.worst, plot=FALSE )
print( all( mcmcResult1$param$beta == mcmcResult2$param$beta ) )
print( all( mcmcResult1$param$sigma == mcmcResult2$param$sigma ) )
print( all( acf.worst1$acf == acf.worst2$acf ) )


## re-run identical mcmc series using data matrices as arguments
Y <- cbind(simdata$Y1, simdata$Y2, simdata$Y3 )
X <- simdata$X1
set.seed(567)
mcmcResult3 <- mcmcNorm( 
   Y, X, starting.values=emResult$param,
   save.worst.series=TRUE, worst.linear.coef=emResult$worst.linear.coef )
acf.worst3 <- acf( mcmcResult3$series.worst, plot=FALSE )
print( all( mcmcResult1$param$beta == mcmcResult3$param$beta ) )
print( all( mcmcResult1$param$sigma == mcmcResult3$param$sigma ) )
print( all( acf.worst1$acf == acf.worst3$acf ) )

# generate five imputations for cholesterol data
data(cholesterol)
emResult <- emNorm(cholesterol)
set.seed(999)
mcmcResult1 <- mcmcNorm(emResult, iter=1000, impute.every=200)
print( round( mcmcResult1$imp.list[[3]], 2) )  # print the third one

## re-run identical series, retaining a one-in-ten subsample
set.seed(999)
mcmcResult2 <- mcmcNorm(emResult, iter=100, multicycle=10, impute.every=20)
print( all( mcmcResult1$imp.list[[3]] == mcmcResult2$imp.list[[3]] ) )

## run MCMC on marijuana data with ridge prior
data(marijuana)
emResult <- emNorm( marijuana,  prior="ridge", prior.df=0.5 )  
set.seed(876)
mcmcResult <- mcmcNorm( emResult )
acf.worst <- acf( mcmcResult$series.worst, plot=FALSE )
print( round( acf.worst$acf, 5 ) )

## run MCMC on flas data with ridge prior
data(flas)
emResult <- emNorm( flas, prior="ridge", prior.df=5 )
set.seed(754)
mcmcResult <- mcmcNorm( emResult )
print( round( mcmcResult$y.imp, 3 ) )
