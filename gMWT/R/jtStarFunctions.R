# Version: 30-11-2012, Daniel Fischer

jtStar <- function(X,g){
  goi <- sort(unique(g))
  diffProbs <- getComb(goi,"pairs",order=T)

  result <- 0
  for(sumRun in 1:dim(diffProbs)[1])
  {
    result <- result + getP.Cnaive(X[g==diffProbs[sumRun,1]],X[g==diffProbs[sumRun,2]])
  }
  result
}

jtStarPTest <- function(X,g,nper){
    permTestValues <- c(rep(0,nper))
    for(i in 1:nper){
	permGroups <- sample(g)
	permTestValues[i] <- jtStar(X,permGroups)
    }
    permTestValues
}