toDirect <- function(adjMat) {
  theNum <- nrow(adjMat)
  for (i in 1:(nrow(adjMat) - 1)) {
    b <- i + 1
    for (j in b:nrow(adjMat)) {
      if (adjMat[i, j] && adjMat[j, i] == 1) {
        adjMat[i, j] <- adjMat[j, i] <- 0
      }
    }
  }
  #so as to have the same matrix format
  return(matrix(adjMat, theNum, theNum))
}


toSkeleton<- function(adjMat) {
  for(i in 1:(nrow(adjMat) - 1)) {
    b <- i + 1
    for (j in b:nrow(adjMat)) {
      if (adjMat[i, j] != adjMat[j, i]) {
        adjMat[i, j] <- adjMat[j, i] <- 1
      }
    }
  }
  #so as to have the same matrix format
  return(matrix(adjMat, nrow(adjMat), nrow(adjMat)))
}


edgeCounter <- function(adjMat) {
  totalEdge <- matrix(0, nrow(adjMat), nrow(adjMat))
  for(i in 1:(nrow(adjMat) - 1)) {
    b <- i + 1
    for(j in b:nrow(adjMat)) {
      if (adjMat[i, j] && adjMat[j, i] == 1) {
        totalEdge[i, j] <- totalEdge[j, i] <- totalEdge[i, j] + 1
      }
    }
  }
  return(totalEdge)
}


causalCounter <- function(adjMat) {
  totalCausal <- matrix(0, nrow(adjMat), nrow(adjMat))
  for (i in 1:nrow(adjMat)) {
    if (any(matrixcalc::matrix.power(adjMat, i) > 0)) {
      totalCausal <- totalCausal +
        #matrixcalc::matrix.power(adjMat, i)
        #so as to have the same matrix format
        matrix(matrixcalc::matrix.power(adjMat, i), nrow(adjMat), nrow(adjMat))
    }
  }

  totalCausal[totalCausal > 0] <- 1

  return(totalCausal)
}
