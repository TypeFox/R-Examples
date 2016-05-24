mutation <- function(modelString, mutRate, numVar, cons, longitudinal){

  if (longitudinal) {
    if (!all(cons == 0)) {
      cons <- cons + (numVar * numVar)
    }
  }

  for (i in 1:length(modelString)) {

    toss <- runif(1, 0, 1)
      # if toss <= mutation rate and i is not constraint
      if (toss <= mutRate && all(cons != i)) {
          modelString[i] <- (modelString[i] - 1)^2
      }
  }

  #if the model is cyclic
  if (ggm::isAcyclic(stringToMatrix1(modelString, numVar,
                                     longitudinal)) != TRUE) {

      theModel <- cycleRepair(stringToMatrix1(modelString,
                                              numVar, longitudinal))
      diag(theModel) <- NA
      intra <- as.vector(theModel)

      if(longitudinal) {
        modelString <- c(modelString[1:(numVar * numVar)],
                         intra[!is.na(intra)])
      } else {
        modelString <- intra[!is.na(intra)]
      }

  }
  return(modelString)
}
