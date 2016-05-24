crossOver <- function(stringModelA, stringModelB, numVar, longitudinal) {
  if (longitudinal) {
    #get intra block from parenta A and B
    #so for longitudinal model, the swap is between intra
    fromA <- stringModelA[(numVar * numVar + 1):length(stringModelA)]
    fromB <- stringModelB[(numVar * numVar + 1):length(stringModelB)]

    #swap to reproduce offsprings
    offSpringA <- c(stringModelA[1:(numVar * numVar)], fromB)
    offSpringB <- c(stringModelB[1:(numVar * numVar)], fromA)

  } else {
    #if cross-sectional data

    # get a half from each parent
    fromA <- stringModelA[1:round(length(stringModelA) / 2)]
    fromB <- stringModelB[1:round(length(stringModelB) / 2)]

    #swap
    offSpringA <- c(fromB, stringModelA[(round(length(stringModelA) /
                                                 2) + 1):length(stringModelA)])
    offSpringB <- c(fromA, stringModelB[(round(length(stringModelB) /
                                                 2) + 1):length(stringModelB)])

    #check whether the children are cyclic, if so then repair
    if(!ggm::isAcyclic(stringToMatrix1(offSpringA, numVar, longitudinal))){

      theModelA <- cycleRepair(stringToMatrix1(offSpringA, numVar, longitudinal))

      diag(theModelA) <- NA
      intraA <- as.vector(theModelA)
      offSpringA <- intraA[!is.na(intraA)]
    }

    if(!ggm::isAcyclic(stringToMatrix1(offSpringB, numVar, longitudinal))){
      theModelB <- cycleRepair(stringToMatrix1(offSpringB, numVar,
                                               longitudinal))

      diag(theModelB) <- NA
      intraB <- as.vector(theModelB)
      offSpringB <- intraB[!is.na(intraB)]
    }
  }
  return(list(offSpringA, offSpringB))
}
