# Test adapted from fields package, under GPL license

library( fields )
options( echo=FALSE)
test.for.zero.flag<- 1

#
##### generate test data
#

genCovMat = function(x, theta, lambda) {
  distanceMatrix<- rdist(x,x)
  Sigma<- exp( -distanceMatrix/theta ) + diag(x=lambda, nrow=nrow(distanceMatrix))
  return(Sigma)
}

#generate observation locations
n=500
x = matrix(runif(2*n), nrow=n)

#generate observations at the locations
trueTheta = .2
trueLambda = .1
Sigma = genCovMat(x, trueTheta, trueLambda)

U = chol(Sigma)
y = t(U)%*%as.vector(rnorm(n))

#
######set MLE computation parameters
#

testThetas = seq(from=trueTheta/2, to=2*trueTheta, length=20)
par.grid=list(theta=testThetas)
guessLambda = trueLambda

#
##### test using distance matrix
#

print("testing using distance matrix")

set.seed(1)
out1 = mKrig.MLE(x, y, lambda=guessLambda, par.grid=par.grid, cov.args= list(Distance="rdist"))
lambda.MLE = out1$lambda.MLE
theta.MLE = out1$cov.args.MLE$theta

#perform mKrig at MLE parameters
out1 = mKrig(x, y, lambda=lambda.MLE, theta=theta.MLE, cov.args= list(Distance="rdist"))
print("finished default case")

set.seed(1)
out2 = mKrig.MLE(x, y, lambda=guessLambda, par.grid=par.grid)
lambda.MLE = out2$lambda.MLE
theta.MLE = out2$cov.args.MLE$theta

#perform mKrig at MLE parameters
out2 = mKrig(x, y, lambda=lambda.MLE, theta=theta.MLE)

print("finished compact distance matrix case")

#
##### test comatibility with other fields functions
#

temp1<- predict( out1)
temp2<- predict( out2)
test.for.zero( temp1, temp2, tag="predict compatibility: rdist with compact versus normal rdist")

#
##### test SE
#

temp1 = predictSE(out1)
temp2 = predictSE(out2)

test.for.zero( temp1, temp2, tag="predictSE compatibility: rdist with compact versus normal rdist")





cat("all done with mKrig.MLE tests", fill=TRUE)
options( echo=TRUE)
