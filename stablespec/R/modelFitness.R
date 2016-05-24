gatherFitness <- function(dataSet, allModelString, sizeSubset,
                         numVar, index, longitudinal, co) {

  result <- matrix(0, nrow(allModelString), 4)

  if (longitudinal) {
    for (i in 1:nrow(allModelString)) {

      theModel <- longiMatrixFill(allModelString[i, ], numVar)
      #convert model into string of characters to be read by sem::specifyEquations
      modelChar <- writeModel(theModel, numVar, longitudinal)
      result[i, ] <- getFitness(dataSet, modelChar, sizeSubset, index, co)
    }

  } else {
    #if cross-sectional
    for (i in 1:nrow(allModelString)) {

      theModel <- stringToMatrix1(allModelString[i, ], numVar, longitudinal)
      #convert model into string of characters to be read by sem::specifyEquations
      modelChar <- writeModel(theModel, numVar, longitudinal)
      result[i, ] <- getFitness(dataSet, modelChar, sizeSubset, index, co)
    }
  }

  return(result)
}


getFitness <- function(dataSet, modelChar, sizeSubset, index, co) {
  # write the model into text file
  model_spec <- suppressMessages(sem::specifyEquations(text = modelChar))

  # a variables for new names of the data set
  theNames <- list()
  for(i in 1:ncol(dataSet)) {
    theNames[[i]] <- paste('var', i, sep="")
  }

  if(co == "covariance") {
    #covariance
    coData <- cov(dataSet)
  } else { #correlation
    coData <- cor(dataSet)
  }

  rownames(coData) <- colnames(coData) <- c(unlist(theNames))
  #compute the SEM
  theFitness <- sem::sem(model_spec, coData, sizeSubset)
  #get the summary
  fitSummary <- summary(theFitness)

  chi <- fitSummary$chisq
  df <- fitSummary$df
  bic <-fitSummary$BIC

  #adhoc handling
  if (df == 0 & chi != 0) {
    chi <- 0
  }
  return(c(chi, df, bic, index))
}
