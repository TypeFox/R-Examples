fastNonDominatedSort <- function (matrixOfFitness){
  # All credit goes to Ching-Shih Tsou who originally wrote this function.
  # We borrow this function and modified it so as to meet our need.
  #
  # Ching-Shih Tsou (2013). nsga2R: Elitist Non-dominated Sorting
  # Genetic Algorithm based on R. R package version 1.0.
  # https://CRAN.R-project.org/package=nsga2R


  popSize <- nrow(matrixOfFitness)
  idxDominators <- vector("list", popSize)
  idxDominatees <- vector("list", popSize)
  for (i in 1:(popSize - 1)) {
    for (j in i:popSize) {
      if (i != j) {
        xi <- matrixOfFitness[i, 1:2]
        xj <- matrixOfFitness[j, 1:2]
        #xi[1] is chisq, xi[2]is df
        if ((xi[1] <= xj[1] && xi[2] > xj[2]) || (xi[1] < xj[1] && xi[2] >= xj[2])) {
          idxDominators[[j]] <- c(idxDominators[[j]],i)
          idxDominatees[[i]] <- c(idxDominatees[[i]],j)
        } else if ((xj[1] <= xi[1] && xj[2] > xi[2]) || (xj[1] < xi[1] && xj[2] >= xi[2])) {
          idxDominators[[i]] <- c(idxDominators[[i]],j)
          idxDominatees[[j]] <- c(idxDominatees[[j]],i)
        }
      }
    }
  }
  noDominators <- lapply(idxDominators, length)
  rnkList <- list()
  rnkList <- c(rnkList, list(which(noDominators == 0)))
  solAssigned <- c()
  solAssigned <- c(solAssigned, length(which(noDominators == 0)))

  while (sum(solAssigned) < popSize) {
    Q <- c()
    noSolInCurrFrnt <- solAssigned[length(solAssigned)]
    for (i in 1:noSolInCurrFrnt) {
      solIdx <- rnkList[[length(rnkList)]][i]
      hisDominatees <- idxDominatees[[solIdx]]
      for (i in hisDominatees) {
        noDominators[[i]] <- noDominators[[i]] - 1
        if (noDominators[[i]] == 0) {
          Q <- c(Q, i)
        }
      }
    }
    rnkList <- c(rnkList, list(sort(Q)))
    solAssigned <- c(solAssigned, length(Q))
  }
  return(rnkList)
}

convertFront <- function(sortedRnk) {
  # Convert list of sorted ranking into matrix
  return(matrix(unlist(sortedRnk), length(unlist(sortedRnk)), 1))
}

sortBasedOnDist <- function(rnk, dist) {
  for(i in 1:length(rnk)) {
    rnk[[i]] <- rnk[[i]][order(-dist[rnk[[i]], ])]
  }
  return(rnk)
}


