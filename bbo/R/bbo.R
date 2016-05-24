bbo.control <- function(pModify = 1, pMutate = 0.3, KEEP = 5, popSize = 20, maxGen = 20, numVar = 2, orderDep = TRUE) {


  if (pModify > 1 | pModify < 0) {
    warning("'pModify' not in [0,1]; set to default value 1\n", immediate. = TRUE)
    pModify <- 1
  }

  if (pMutate < 0 | pMutate > 1) {
    warning("'pMutate' not in [0,1]; set to default value 0.3\n", immediate. = TRUE)
    pMutate <- 0.3
  }

  if (popSize <= 0 | popSize < KEEP) {
    warning("'popSize' must be > 0 and/or greater than 'KEEP'; 'popSize' set to default value 20\n", immediate. = TRUE)
    #popSize <- 20
    stop("Halted!")
  }

   if (KEEP < 0 | KEEP > popSize) {
    warning("'KEEP' cannot be negative or greater than 'popSize'; 'KEEP' default value is 5\n", immediate. = TRUE)
    stop("incorrect value passed!")
  }

  if (maxGen <= 0) {
    warning("'maxGen' must be > 0, set to default value 20\n", immediate. = TRUE)
    maxGen <- 20
  }


  if (numVar <= 0) {
    warning("'numVar' cannot be negative or O; numVar default value is 2\n", immediate. = TRUE)
    stop("incorrect value passed!")
  }

  if (missing(orderDep)){
    #warning("'orderDep' is set to TRUE\n", immediate. = TRUE)
    orderDep <- TRUE
  }


  list(pModify = pModify, pMutate = pMutate, KEEP = KEEP, popSize = popSize, maxGen = maxGen, numVar = numVar, orderDep = orderDep)

}


bbo <- function(fn, lower, upper, DisplayFlag = TRUE, RandSeed, control = bbo.control(), ... ) {
  env <- new.env()
  ctrl <- do.call(bbo.control, as.list(control))
  selfSet <- FALSE
  if(missing(RandSeed)){
    RandSeed <- 1024
    selfSet <- TRUE
  }

  lower <- rep(lower, ctrl$numVar)
  upper <- rep(upper, ctrl$numVar)

  if (length(lower) != length(upper))
    stop("'lower' and 'upper' are not of same length")

  if (!is.vector(lower))
    lower <- as.vector(lower)

  if (!is.vector(upper))
    upper <- as.vector(upper)

  if (any(lower > upper))
    stop("'lower' > 'upper'")

  if (any(lower == "Inf"))
    warning("A component of 'lower' assigned 'Inf'. May imply 'NaN' results", immediate. = TRUE)

  if (any(lower == "-Inf"))
    warning("A component of 'lower' assigned '-Inf'. May imply 'NaN' results", immediate. = TRUE)

  if (any(upper == "Inf"))
    warning("A component of 'upper' assigned 'Inf'. May imply 'NaN' results", immediate. = TRUE)

  if (any(upper == "-Inf"))
    warning("A component of 'upper' assigned '-Inf'. May imply 'NaN' results", immediate. = TRUE)


  #if(length(lower) != length(upper)){
  #  warning("length of vectors 'lower' and 'upper' have to be equal.", immediate. = TRUE)
  #  stop("Halted!")
  #}
  #if(length(lower) > 1 & length(lower) != ctrl$popSize){
  #  warning("length of vectors 'lower' and 'upper' must be equal to the initial population size")
  #  stop("Halted!")
  #}

  #if(length(lower) == 1){
  #  lower <- rep(lower, ctrl$numVar)
  #  upper <- rep(upper, ctrl$numVar)
  #}
  
  if(!missing(RandSeed)){
    if(!selfSet){
      set.seed(RandSeed)
    }
  }

  ## FeasibleFunction = #todo (Checks whether each candidate solution is a feasible solution; if not it bounds the individual parameter values by the minimum and maximum)

  ## Supporting functions::
  #TO DO# ClearDups : makes sure that the population does not have duplicates


  ## InitPopulation function code
  ## Next version:  we may allow the user to supply his own set of initial population, if that seems feasible considering all possible constraints

  members <- matrix(runif(ctrl$popSize * ctrl$numVar, lower, upper), nrow = ctrl$popSize, ncol = ctrl$numVar, byrow = TRUE)

  ## Pass the additional arguments to the objective function
  costs <- apply(members, 1, fn, ...)

  ## 'members' is a matrix, 'costs' is a vector

  population <- list(members = members, costs = costs)

  eliteSolutions <- matrix(rep(-Inf,ctrl$KEEP * ctrl$numVar), nrow = ctrl$KEEP, ncol = ctrl$numVar)
  eliteCosts <- rep(Inf, ctrl$KEEP)

  MinCostEachGen <- c()
  AvgCostEachGen <- c()
  BestMemberEachGen <- c()

  ## Begin the optimization loop
  for( genIndex in 1:ctrl$maxGen ){
    
    ## Save the best habitats in a temporary array/vector
  
    eliteSolutions <- as.matrix(population$members[1:ctrl$KEEP, ])
    eliteCosts <- population$costs[1:ctrl$KEEP]
    #elite <- list(eliteSolutions, eliteCosts)

    ## Compute immigration rate and emigration rate for each species count.
    ## lambda(i) is the immigration rate for habitat i.
    ## mu(i) is the emigration rate for habitat i.
    
    lambda <- rep(Inf, ctrl$popSize)
    mu <- rep(Inf, ctrl$popSize)

    ## mu(i) is the extinction rate for individual i.
    ## This routine assumes the population is sorted from most fit to least fit.
    ## getLambdaMu() is a separate function in Matlab impelementation, while here we do it in-place
    
    P <- nrow(population$members);
    for( i in 1:P ){
      if(population$costs[i] < Inf){
        mu[i] <- (P - i) / P
      }else{
        mu[i] <- 0
      }
    lambda[i] <- 1 - mu[i]
    }
    ## getLambdaMu done


    lambdaMin <- min(lambda)
    lambdaMax <- max(lambda)

    Island <- matrix(-Inf, nrow = nrow(population$members), ncol = ncol(population$members))

    for( k in 1:nrow(population$members) ){
        if(runif(1) > ctrl$pModify){
            ## do-nothing
            next
        }
	
        ## Normalize the immigration rate.
        lambdaScale <- (lambda[k] - lambdaMin) / (lambdaMax - lambdaMin)

        ## Probabilistically input new information into habitat i
        for( j in 1:ctrl$numVar ){
            if(runif(1) < lambdaScale){
                ## Pick a habitat from which to obtain a feature
                RandomNum <- runif(1) * sum(mu)
                Select <- mu[1]
                SelectIndex <- 1
                while((RandomNum > Select) & (SelectIndex < ctrl$popSize)){
                    SelectIndex <- SelectIndex + 1
                    Select <-Select + mu[SelectIndex]
                }## while-ends
                Island[k,j] <- population$members[SelectIndex,j]
	    }
            else{
                Island[k,j] <- population$members[k,j]
            }## if-else-ends
        }## second-for-ends
    }## first-for-ends

    ## Mutation ##
    population <- popSort(population)
    for( k in 1:nrow(population$members) ){
        for( parnum in 1:ctrl$numVar ){
            if( ctrl$pMutate > runif(1) ){
                Island[k,parnum] <- runif(1, lower[parnum], upper[parnum])
            }
        }
    }

    ## Replace the habitats with their new versions.
    for( k in 1:nrow(population$members)){
        population$members[k, ] <- Island[k, ]
    }

    ## Make sure each individual is legal.
    population <- feasibilityCheck(population, lower, upper)
    
    ## Calculate cost
    population$costs <- apply(population$members, 1, fn)

    ## Sort from best to worst
    population <- popSort(population)
    
    ## Replace the worst with the previous generation's elites.
    n <- nrow(population$members)

    if(ctrl$KEEP > 0){
      for( k in 1:ctrl$KEEP ){
          population$members[n-k+1, ] <- eliteSolutions[k, ]
          population$costs[n-k+1] <- eliteCosts[k]
      }
    }

    ## Make sure the population does not have duplicates. 
    ## Population <- clearDups(population, upper, lower, fn)
    
    ## Sort from best to worst
    population <- popSort(population)
    
    ## Compute the average cost
    ## nLegal: return value from computeAveCost is not used anywhere, so we dont collect it, for now.
    AverageCost <- computeAveCost(population)

    ## Display info to screen
    MinCostEachGen <- c(MinCostEachGen, population$cost[1])
    AvgCostEachGen <- c(AvgCostEachGen, AverageCost)
    currBestPop <- population$members[1, ]
    BestMemberEachGen <- rbind(BestMemberEachGen, currBestPop)
    rownames(BestMemberEachGen) <- NULL
    if( DisplayFlag ){
        cat('The best and mean of Generation # ', genIndex, ' are ', tail(MinCostEachGen, 1), ' and ', tail(AvgCostEachGen, 1), '\n')
    }


  }## main-optimization-for-loop-upto-maxgen-ends

  bestMember <- population$members[1, ] 
  bestValue <-  tail(MinCostEachGen, 1)
  ## conclude function can be called here
  ## FINAL RETURN
  minCost <- list(bestMember = bestMember, bestValue = bestValue)
  bestHabitat <- list(itersBestValue = MinCostEachGen, itersBestMember = BestMemberEachGen, itersAvg = AvgCostEachGen)

  return(list(prop = ctrl, minCost = minCost, bestHabitat = bestHabitat))

}


feasibilityCheck <- function(population, lower, upper){

for( i in 1:nrow(population$members)){
    for( k in 1:ncol(population$members)){
        population$members[i, k] <- max(population$members[i, k], lower[k])
        population$members[i, k] <- min(population$members[i, k], upper[k])
    }
}

return ( population )
}


computeAveCost <- function(population){

## Compute the average cost of all legal individuals in the population.
## OUTPUTS: avgCost = average cost
##          nLegal = number of legal individuals in population
## Save valid population member fitnesses in temporary array

Cost <- c(); nLegal <- 0;

for( i in 1:length(population$costs)){
    if( population$costs[i] < Inf ){
        Cost <- c(Cost, population$costs[i])
        nLegal <- nLegal + 1
    }
}

## Compute average cost.
avgCost <- mean(Cost)
return ( avgCost )
}

popSort <- function(population){

##Sort the population members from best to worst

population$members <- as.matrix(population$members[order(population$costs), ])
population$costs <- population$costs[order(population$costs)]
return ( population )
}


## This function is kept on hold!!
##clearDups <- function(population, upper, lower, fn){
##return ( population )
##}
##
