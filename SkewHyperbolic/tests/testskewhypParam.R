require(SkewHyperbolic)
## source("../R/skewhypCalcRange.R")
## source("../R/dskewhyp.R")
## source("../data/skewhypParam.R")

### Testing the accuracy of skewhypMean
data(skewhypParam)
for (i in 1:nrow(skewhypSmallParam)) {
  param <- skewhypSmallParam[i, ]
  x <- rskewhyp(1000, param = param)
  sampleMean <- mean(x)
  distMean <- skewhypMean(param = param)
  difference <- abs(sampleMean - distMean)
  print(difference)
}
