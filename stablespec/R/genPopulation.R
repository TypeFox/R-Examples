genPopulation <- function(nPop, numVar, longitudinal, consVector) {
  modelPopulation <- NULL

  for (i in 1:nPop) {
    # generate random string for intra-slice
    intraString <- round(runif(numVar * (numVar - 1)))
    #subtitution of   stringModel = round(runif(numVar*(numVar-1)))

    # satisfy the constraint
    if (any(intraString[consVector] != 0)) {
      intraString[consVector] <- 0
    }


    if (longitudinal) {
      #generate random string for inter-slice
      interString <- round(runif(numVar * numVar))
      stringModel <- c(interString, intraString)

      #check whether this new string model is cyclic
      tmpModel <- stringToMatrix1(stringModel, numVar, longitudinal)
      if (!ggm::isAcyclic(tmpModel)) {
        stringModel <- repairCyclicModel(stringModel, numVar, longitudinal)
      }

      #stringModel <- repairCyclicModel(stringModel, numVar, longitudinal)
      theModel <- matrix(stringModel, 1,
                         (numVar * (numVar - 1)) + (numVar * numVar))
      #the length of the string is intra + inter-slice
    } else {

      #if cross-sectional
      stringModel <- intraString

      #check whether this new string model is cyclic
      tmpModel <- stringToMatrix1(stringModel, numVar, longitudinal)
      if (!ggm::isAcyclic(tmpModel)) {
        stringModel <- repairCyclicModel(stringModel, numVar, longitudinal)
      }

      #stringModel <- repairCyclicModel(intraString, numVar, longitudinal)
      theModel <- matrix(stringModel, 1, numVar * (numVar - 1))
    }
    modelPopulation <- rbind(modelPopulation, theModel)
  }
  return(modelPopulation)
}
