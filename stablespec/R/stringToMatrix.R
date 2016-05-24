stringToMatrix1 <- function(stringModel, numVar, longitudinal) {
  if (longitudinal) {

    #get only the intra string
    stringModel_intra <- stringModel[(numVar * numVar + 1):length(stringModel)]

  } else { #if cross-sectional

    stringModel_intra <- stringModel
  }

    theModel <- matrix(0, numVar, numVar)

    for (i in 1:numVar) {
      for (j in 1:numVar) {
        if (i != j ) {
          theModel[j, i] <- head(stringModel_intra, n=1)
          stringModel_intra <- stringModel_intra[-1]
        }
      }
    }

  return(theModel)
}



longiMatrixFill <- function(stringModel, numVar) {
  # for filling in longitudinal matrix, basically the same concept
  # as above function but this function is mainly built for filling,
  # while above function can be also used to check
  # the intra part of longitudinal model, eg cyclic check; so faster.

  theModel <- matrix(0, (2 * numVar), (2 * numVar))
  stringModel_intra <- stringModel[(numVar * numVar + 1):length(stringModel)]
  stringModel_inter <- stringModel[1:(numVar * numVar)]


  #fill in intra matrix
  the_loop1 <- (1:numVar) + numVar

  for (i in the_loop1) {
    for (j in the_loop1) {
      if (i != j) {
        theModel[j, i] <- head(stringModel_intra, n=1)
        stringModel_intra <- stringModel_intra[-1]
      }
    }
  }

  #fill in inter matrix
  the_loop1 <- 1:numVar
  the_loop2 <- (1:numVar) + numVar

  for (i in the_loop2) {
    for (j in the_loop1) {
      theModel[j, i] <- head(stringModel_inter, n=1)
      stringModel_inter <- stringModel_inter[-1]
    }
  }

  return(theModel)
}
