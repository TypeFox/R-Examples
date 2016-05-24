library(HyperbolicDist)
if (FALSE){
  fileList <- list.files("../R/")
  fileList <- fileList[-grep(".+~",fileList)]
  fileNames <- paste("../R/", fileList, sep = "")
  lapply(fileNames,source)
}
sampSize <- 100000
Theta <- c(-0.5,5,2.5)
dataVector <- rgig(sampSize, Theta)
### Test gigMean
(m1 <- gigMean(Theta))
momIntegrated("gig", order = 1, param = Theta, about = 0)
mean(dataVector)

### Test gigVar
gigVar(Theta)
gigMom(2, Theta, about = m1)
momIntegrated("gig", order = 2, param = Theta, about = m1)
var(dataVector)

### Test gigSkew
gigSkew(Theta)
skewness(dataVector)

### Test gigKurt
gigKurt(Theta)
kurtosis(dataVector)


