library(HyperbolicDist)
if (FALSE){
  fileList <- list.files("../R/")
  fileList <- fileList[-grep(".+~",fileList)]
  fileNames <- paste("../R/", fileList, sep = "")
  lapply(fileNames,source)
}

Theta <- c(2,2,1,2,1)
mu <- Theta[5]
### Test ghypMom: mu moments
(m1 <- ghypMean(Theta))
m1 - mu
ghypMom(1, Theta, momType = "mu")
momIntegrated("ghyp", order = 1, param = Theta, about = mu)
ghypMom(2, Theta, momType = "mu")
momIntegrated("ghyp", order = 2, param = Theta, about = mu)
ghypMom(10, Theta, momType = "mu")
momIntegrated("ghyp", order = 10, param = Theta, about = mu)

### Test ghypMom: raw moments
ghypMean(Theta)
ghypMom(1, Theta, momType = "raw")
momIntegrated("ghyp", order = 1, param = Theta, about = 0)
ghypMom(2, Theta, momType = "raw")
momIntegrated("ghyp", order = 2, param = Theta, about = 0)
ghypMom(10, Theta, momType = "raw")
momIntegrated("ghyp", order = 10, param = Theta, about = 0)

### Test ghypMom: central moments
ghypMom(1, Theta, momType = "central")
momIntegrated("ghyp", order = 1, param = Theta, about = m1)
ghypVar(Theta)
ghypMom(2, Theta, momType = "central")
momIntegrated("ghyp", order = 2, param = Theta, about = m1)
ghypMom(10, Theta, momType = "central")
momIntegrated("ghyp", order = 10, param = Theta, about = m1)


