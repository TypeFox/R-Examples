library(HyperbolicDist)
if (FALSE){
  fileList <- list.files("../R/")
  fileList <- fileList[-grep(".+~",fileList)]
  fileNames <- paste("../R/", fileList, sep = "")
  lapply(fileNames,source)
}
sampSize <- 100000
Theta <- c(2,2,1,2,1)
mu <- Theta[5]
dataVector <- rghyp(sampSize, Theta)
### Test ghypMean
(m1 <- ghypMean(Theta))
momIntegrated("ghyp", order = 1, param = Theta, about = 0)

### Test ghypVar
ghypVar(Theta)
ghypMom(2, Theta, momType = "central")
momIntegrated("ghyp", order = 2, param = Theta, about = m1)
var(dataVector)

### Test ghypSkew
ghypSkew(Theta)
skewness(dataVector)

### Test ghypKurt
ghypKurt(Theta)
kurtosis(dataVector)


