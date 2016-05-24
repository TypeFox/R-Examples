structureStab <- function(listOfFronts, stringSize, numVar,
                          longitudinal, consMatrix) {

  #convert the consMatrix
  consMatrix <- cons4Stab(consMatrix, numVar, longitudinal)

  # Initialization of lists
  listOfCausal <- listOfCausal_l1 <- listOfEdge <- list()

  # To get only a unique model (of a front) in each model complexity
  for (i in 1:length(listOfFronts))
  {
    tempMatA <- listOfFronts[[i]]
    listOfFronts[[i]] <- tempMatA[match(unique(tempMatA[, stringSize + 2]),
                                        tempMatA[, stringSize + 2]), ]
  }

  # Convert list of optimal models into a matrix
  matOfFronts <- do.call(rbind, listOfFronts)

  # To list optimal models based on model complexity
  groupedFronts <- groupFrontComp(matOfFronts, stringSize)


  if (longitudinal) {

    a <- min(matOfFronts[, stringSize + 2]) + 1
    b <- max(matOfFronts[, stringSize + 2]) + 1

    for (i in a:b) {

      listOfCausal[[i]] <- listOfCausal_l1[[i]] <-
        listOfEdge[[i]] <- matrix(0, 2 * numVar, 2 * numVar)

      if (any(groupedFronts[[i]] > 0)) {

        matModels <- groupedFronts[[i]]

        for (j in 1:nrow(matModels)) {

          tempMatB <- longiMatrixFill(matModels[j, ], numVar)
          theGraph <- as(tempMatB, Class="graphNEL")

          #convert to CPDAG
          theCPDAG <- dag2CpdagCons(theGraph, consMatrix)

          # extract the compelled and sckeleton
          adjMatCausal <- toDirect(as(theCPDAG, "matrix"))
          adjMatSkel <- toSkeleton(as(theCPDAG, "matrix"))

          #count the occurences
          listOfCausal[[i]] <- listOfCausal[[i]] + causalCounter(adjMatCausal)
          listOfCausal_l1[[i]] <- listOfCausal_l1[[i]] + adjMatCausal
          listOfEdge[[i]] <- listOfEdge[[i]] + edgeCounter(adjMatSkel)

        }
        # compute the average
        listOfCausal[[i]] <- round(listOfCausal[[i]] / nrow(matModels), 2)
        listOfCausal_l1[[i]] <- round(listOfCausal_l1[[i]] / nrow(matModels), 2)
        listOfEdge[[i]] <- round(listOfEdge[[i]] / nrow(matModels), 2)
      }
    }
  } else { #if cross-sectional
    b <- max(matOfFronts[, stringSize + 2]) + 1

    for (i in 1:b) {
      listOfCausal[[i]] <- listOfCausal_l1[[i]] <-
        listOfEdge[[i]] <- matrix(0, numVar, numVar)

      if (any(groupedFronts[[i]] > 0)) {

        matModels <- groupedFronts[[i]]

        for (j in 1:nrow(matModels)) {

          tempMatB <- stringToMatrix1(matModels[j, ], numVar, longitudinal)
          theGraph <- as(tempMatB, Class="graphNEL")

          #convert to CPDAG
          theCPDAG <- dag2CpdagCons(theGraph, consMatrix)

          # extract the compelled and sckeleton
          adjMatCausal <- toDirect(as(theCPDAG, "matrix"))
          adjMatSkel <- toSkeleton(as(theCPDAG, "matrix"))

          #count the occurences
          listOfCausal[[i]] <- listOfCausal[[i]] + causalCounter(adjMatCausal)
          listOfCausal_l1[[i]] <- listOfCausal_l1[[i]] + adjMatCausal
          listOfEdge[[i]] <- listOfEdge[[i]] + edgeCounter(adjMatSkel)
        }

        # compute the average
        listOfCausal[[i]] <- round(listOfCausal[[i]] / nrow(matModels), 2)
        listOfCausal_l1[[i]] <- round(listOfCausal_l1[[i]] / nrow(matModels), 2)
        listOfEdge[[i]] <- round(listOfEdge[[i]] / nrow(matModels), 2)
      }
    }
  }

  # to handle adhoc situation
  listOfCausal[sapply(listOfCausal, is.null)] <-
    listOfCausal_l1[sapply(listOfCausal_l1, is.null)] <-
    listOfEdge[sapply(listOfEdge, is.null)] <- NULL

  return(list(causalStab=listOfCausal, causalStab_l1=listOfCausal_l1,
              edgeStab=listOfEdge))
}
