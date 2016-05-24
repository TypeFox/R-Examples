convertCons <- function(consMatrix, numVar) {
  # create temporary matrix
  tempMat <- matrix(NA, numVar, numVar)

  # initial index
  ind <- 1

  #fill in the matrix by column
  for (i in 1:numVar) {
    for (j in 1:numVar) {
      if (i != j) {
        tempMat[j, i] <- ind
        ind <- ind + 1
      }
    }
  }

conString <- NULL

  for (i in 1:nrow(consMatrix)) {
    conString <- c(conString, tempMat[consMatrix[i, 1], consMatrix[i, 2]])
  }

  return(sort(conString))
}

#function to convert/extend the constraint matrix for stability selection
cons4Stab <- function(consMatrix, numVar, longitudinal) {

  if (longitudinal) {

    #swap column, for intra, and add with NumVar to get correct indices
    #in longitudinal model
    cons_intra <- consMatrix[, c(2, 1)] + numVar

    # constraint for inter
    cons_inter <- NULL
    for (i in 1:numVar) {
      cons_inter <- rbind(cons_inter,
                          matrix(c(c(1:numVar), rep(numVar + i, numVar)),
                                 numVar, 2))
    }
    return(rbind(cons_intra, cons_inter))

  } else { #if cross-sectional

    #swap the column, e.g., 1 becomes 2 and vice versa
    return(consMatrix[, c(2, 1)])

    }
}
