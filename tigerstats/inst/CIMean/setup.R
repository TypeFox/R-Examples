## Set up the underlying populations for the app


########################################
#  Utility Functions for Generating Pareto Values
########################################

rpareto <- function(n,alpha,theta) {#random values for Pareto(alpha,theta) distribution
  theta*((1-runif(n))^(-1/alpha)-1)
}

dpareto <- function(x,alpha,theta) {  #pdf for Pareto(alpha,theta) distribution
  alpha*theta^alpha/(x+theta)^(alpha+1)
}


#############################################
# Generate the populations
############################################
muNorm <- 70
sigmaNorm <- 5
shapeGamma <- 2
scaleGamma <- 50

# for pareto:
alphaPareto <- 5
thetaPareto <- 100
tailProb <- 0.02  #want to find a Value at risk of 1 - this
valRisk <- thetaPareto*(tailProb^(-.5)-1)

# for pop with group of outliers
propOutliers <- 0.10
meanOutliers <- 200
sdOutliers <- 5
meanRegulars <- 50
sdRegulars <- 5

routlier <- function(n) {
  propNormals <- 1- propOutliers
  whichHump <- rbinom(n,size=1,prob=propNormals)
  outlierSamp <- ifelse(whichHump,rnorm(n,mean=meanRegulars,sd=sdRegulars),
                        rnorm(n,mean=meanOutliers,sd=sdOutliers))
  outlierSamp
}

######################################
# Make population densities
#####################################
xNorm <- seq(muNorm-5*sigmaNorm,muNorm+5*sigmaNorm,length.out=600)
yNorm <- dnorm(xNorm,mean=muNorm,sd=sigmaNorm)
normalDen <- list(x=xNorm,y=yNorm)

xSkew <- seq(0,shapeGamma*scaleGamma+7.5*sqrt(shapeGamma)*scaleGamma,
             length.out=600)
ySkew <- dgamma(xSkew,shape=shapeGamma,scale=scaleGamma)
skewDen <- list(x=xSkew,y=ySkew)

xSuperSkew <- seq(0,valRisk,length.out=600)
ySuperSkew <- dpareto(xSuperSkew,alpha=alphaPareto,theta=thetaPareto)
superSkewDen <- list(x=xSuperSkew,y=ySuperSkew)

xOut <- seq(0,meanOutliers+5*sdOutliers,length.out=600)
yOut <- (1-propOutliers)*dnorm(xOut,mean=meanRegulars,sd=sdRegulars)+propOutliers*dnorm(xOut,mean=meanOutliers,sd=sdOutliers)
outlierDen <- list(x=xOut,y=yOut)

#######################################
# Get the population means
######################################

normalMean <- muNorm
skewMean <- shapeGamma*scaleGamma
superSkewMean <- thetaPareto/(alphaPareto - 1)
outlierMean <- (1-propOutliers)*meanRegulars+propOutliers*meanOutliers