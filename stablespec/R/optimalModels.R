optimalModels <- function(theData, nSubset, iteration, nPop,
                          mutRate, crossRate, longitudinal,
                          numTime, seed, co, consMatrix) {

  #masterList to store the pareto fronts
  #allSubsets to store the subsets
  masterList <- allSubsets <- list()


  if (longitudinal){
    #get the number of instances of one time slice
    numInstances <- nrow(theData) / (numTime - 1)

    #compute the size of Subset
    sizeSubset <- floor(numInstances / 2) * (numTime - 1)

    #compute sizeSubsetGetData, how many instances drawn from each time slice
    sizeSubsetGetData <- floor(numInstances / 2)

    #compute the number of variables,
    #asummed data already in reshaped network t_0..t_i
    numVar <- ncol(theData) / 2

    #compute the size of the a longitudinal model string
    stringSize <- (numVar * numVar + (numVar * (numVar - 1))) #inter + intra

  } else {
    #for variable sizeSubset
    sizeSubset <- floor(nrow(theData) / 2)

    #for variable numVar, the number of variables in the data
    numVar <- ncol(theData)

    #the size of a model string
    stringSize <- (numVar * (numVar - 1))
  }

  #form constraint
  if (is.null(consMatrix)) {
    constraint1 <- 0
  } else {
    constraint1 <- convertCons(consMatrix, numVar)
  }


  # the main loop / outer loop -----------------------------------------

  for (l in 1:nSubset) {

    #set seed
    set.seed(seed[l])

    #draw a subset from the data
    if (longitudinal) {

      newData <- getDataLongi(theData, numTime, sizeSubsetGetData)

    } else { #if cross-sectional data
      newData <- getDataCross(theData, sizeSubset)
    }

    #store into the list
    allSubsets[[l]] <- newData


    #make an initial population
    initialModels <- initialPopulation(numVar, stringSize,
                                       longitudinal, consMatrix)
    initialModels <- rbind(initialModels,
                           genPopulation(nPop-nrow(initialModels), numVar,
                                         longitudinal, constraint1))


    #perform initial computation
  initialResult <- gatherFitness(newData, initialModels, sizeSubset,
                                 numVar, l, longitudinal, co)

  # FIRST GENERATION----------------------------------------------------

  # sorted first generation, here is P0
  preP0 <- fastNonDominatedSort(initialResult)

  #to get the matrix of indices from Front1..n
  indexP0 <- convertFront(preP0)
  P0 <- initialModels[indexP0, ]


  # initial Selection
  currentPop <- P0

  # to get the first front's fitness
  front1Fitness <- initialResult[preP0[[1]], ]

  #to get the first front's models
  front1Model <- initialModels[preP0[[1]], ]


  # the loop of NSGA --------------------------------------------

  for (j in 1:iteration) {
    # Crossover
    i <- 1
    while (i < nPop) {
      toss <- runif(1,0,1)
      if (toss <= crossRate) {
        #generate crossovered offspring, by taking two consecutive models
        offsprings <- crossOver(currentPop[i, ],
                                currentPop[i+1, ], numVar, longitudinal)
        currentPop[i, ] <- offsprings[[1]]
        currentPop[i+1,] <- offsprings[[2]]
      }
      #to take the next two consecutive models
      i <- i + 2
    }

    # MUTATION -----------------------------
    for (i in 1:nPop) {
      currentPop[i,] <- mutation(currentPop[i, ], mutRate,
                                 numVar, constraint1, longitudinal)
    }

    allR0 <- rbind(P0, currentPop)
    preR0 <- gatherFitness(newData, allR0, sizeSubset, numVar,
                           l, longitudinal, co)

    ranking <- fastNonDominatedSort(preR0)
    rnkIndex <- integer(nPop)
    i <- 1
    while (i <= length(ranking)) {
      rnkIndex[ranking[[i]]] <- i
      i <- i + 1
    }


    #to get the range of each objective
    objRange <- apply(preR0[,1:2], 2, max) - apply(preR0[,1:2], 2, min)


    #to assign the crowded distance for each member in each front
    assignedDist <- crowdingDistance(preR0, ranking, objRange, numVar)


    #to sort the members in each front in R0 based on the crowded distance
    sortedDist <- sortBasedOnDist(ranking, assignedDist)

    #order and get R0
    indexR0 <- convertFront(sortedDist)
    R0 <- allR0[indexR0, ]

    #new P0
    P0 <- R0[1:nPop, ]


    #to get the first front fitness
    #started from j+1, because at j,
    #there is initial result before the looping
    if (length(sortedDist[[1]]) < nPop) {
      lengthInd <- length(sortedDist[[1]])
      takenInd <- sortedDist[[1]]
    } else {
      lengthInd <- nPop
      takenInd <- sortedDist[[1]][1:nPop]
    }


    front1Fitness <- rbind(front1Fitness, preR0[takenInd, ])

    #to get the first front models
    front1Model <- rbind(front1Model, allR0[takenInd, ])

    # New currentPop
    currentPop <- cbind(P0, rnkIndex[indexR0[1:nPop, ]],
                        assignedDist[indexR0[1:nPop, ], ])
    currentPop <- nsga2R::tournamentSelection(currentPop, nPop, 2)
    currentPop <- currentPop[, 1:stringSize]
  }

  toFrontOneIndex <- fastNonDominatedSort(front1Fitness)

  #to get the model of 1st front of all FrontOne
  firstFront1Model <- front1Model[toFrontOneIndex[[1]], ]

  #to get the fitness of 1st front of FrontOne
  firstFront1Fitness <- front1Fitness[toFrontOneIndex[[1]], ]

  firstFront1Merged <- cbind(firstFront1Model, firstFront1Fitness)

  uniqueFirstFront <- unique(firstFront1Merged)

  masterList[[l]] <- uniqueFirstFront
  }
  return(list(listOfFronts=masterList, num_var=numVar, string_size=stringSize))
}

