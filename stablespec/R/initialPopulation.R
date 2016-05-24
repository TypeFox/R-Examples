initialPopulation <- function(numVar, stringSize, longitudinal, consMatrix) {

  # how many models representing the whole range of model complexity
  numModels <- (numVar * (numVar + 1 ) / 2) - numVar

  # whether this condition or not; both lower and upper combined
  # thus, that is possibly of cyclic model, need to call repairCycleModel()
  both <- FALSE

  interString <- rep(0, numVar * numVar)

  # matrix of models
  # also will be a frist model which all zero represdenting model with
  # no connection
  allString <- matrix(0, 1, stringSize)

  #first check whether the constraint is at lower or upper diagonal
  lowerMat <- which(lower.tri(matrix(0, numVar, numVar)), arr.ind = TRUE)
  upperMat <- which(upper.tri(matrix(0, numVar, numVar)), arr.ind = TRUE)


  #if constraints belong to both upper and lower diagonal
  if (any(duplicated(rbind(upperMat, consMatrix))) &
      any(duplicated(rbind(lowerMat, consMatrix)))) {

    the_index <- rbind(lowerMat, upperMat)

    #exclude those which matches the constraints
    for (i in 1:nrow(consMatrix)) {
      the_index <-  the_index[-which(the_index[, 1] ==
                                       consMatrix[i, 1] &
                                       the_index[, 2] == consMatrix[i, 2]), ]
    }

    #take only numMOdels of them
    the_index <- the_index[1:numModels, ]

    both <- TRUE

  } else if (any(duplicated(rbind(upperMat, consMatrix)))) {

    #if any constraints belong to upper diagonal matrix
    the_index <- lowerMat

  } else {

    #if any constraints belong to upper lower matrix
    the_index <- upperMat
  }

  # for each complexity, generate a model
  for (i in 1:numModels) {
    model <- matrix(0, numVar, numVar)
    for (j in 1:i) {
      model[the_index[j, 1], the_index[j, 2]] <- 1
    }

    #to convert back from matrix to a binary string and bind into allString
    diag(model) <- NA
    intraString <- as.vector(model)
    intraString <- intraString[!is.na(intraString)]

    if (longitudinal) {
      stringModel <- c(interString, intraString)

      if (both) {
        #then repair, as it is possibly cyclic
        stringModel <- repairCyclicModel(stringModel, numVar, longitudinal)
        allString <- rbind(allString, matrix(stringModel, 1, stringSize))

      } else {
        allString <- rbind(allString, matrix(stringModel, 1, stringSize))
      }

    } else { # if cross-sectional
      if (both) {
        #then repair, as it is possibly cyclic
        stringModel <- repairCyclicModel(intraString, numVar, longitudinal)
        allString <- rbind(allString, matrix(stringModel, 1, stringSize))

      } else {
        allString <- rbind(allString, matrix(intraString, 1, stringSize))
      }
    }
  }

  # if longitudinal data, combine the following random
  # interString with the last intraString
  # (represent the most possible complex model
  # in intra-slice relationships)
  if (longitudinal) {
    for(i in 1:(numVar * numVar)) {

      interString <- rep(0, numVar * numVar)

      #take i random index
      the_index <- sample(1:(numVar * numVar), i)
      interString[the_index] <- 1
      stringModel <- c(interString, intraString)
      allString <- rbind(allString, matrix(stringModel, 1, stringSize))
    }
  }
  return(allString)
}
