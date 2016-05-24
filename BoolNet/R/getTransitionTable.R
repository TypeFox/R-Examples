# Retrieves the transition table of an AttractorInfo structure
# in a readable way
getTransitionTable <- function(attractorInfo)
{
  stopifnot(inherits(attractorInfo,"AttractorInfo") || inherits(attractorInfo,"SymbolicSimulation"))
  
  if (inherits(attractorInfo,"SymbolicSimulation"))
  {
    if (is.null(attractorInfo$graph))
      stop(paste("This SymbolicSimulation structure does not contain transition table information.",
           "Please re-run simulateSymbolicModel() with returnGraph=TRUE!"))
    return(attractorInfo$graph)           
  }
  else
  if (is.null(attractorInfo$stateInfo$table))
    stop(paste("This AttractorInfo structure does not contain transition table information.",
           "Please re-run getAttractors() with a synchronous search and returnTable=TRUE!"))

  return(.internal.getTransitionTable(attractorInfo$stateInfo))
}

# Internally used to extract the information from the <stateInfo> subcomponent
# of AttractorInfo
.internal.getTransitionTable <- function(stateInfo)
{
  fixedGenes <- which(stateInfo$fixedGenes != -1)
  nonFixedGenes <- which(stateInfo$fixedGenes == -1)
  
  if (!is.null(stateInfo$initialStates))
  {
    initialStates <- t(apply(stateInfo$initialStates,2,dec2bin,length(stateInfo$genes)))
  }
  else
  {
    initialStates <- matrix(ncol=length(stateInfo$genes),nrow=2^length(nonFixedGenes))
  
    # encode the initial states
    # first, encode the changing part by calculating all combinations
    temp <- allcombn(2,length(nonFixedGenes)) - 1
    initialStates[,nonFixedGenes] <- temp[,ncol(temp):1]
  
    if (length(fixedGenes) > 0)
    # if there are fixed genes, encode their values
      initialStates[,fixedGenes] <- sapply(stateInfo$fixedGenes[fixedGenes],
                   function(value)rep(value,2^length(nonFixedGenes)))
  }
                   
  # build return structure of class TransitionTable
  res <- data.frame(initialStates,
        t(apply(stateInfo$table,2,dec2bin,length(stateInfo$genes))),
        stateInfo$attractorAssignment,stateInfo$stepsToAttractor)
  colnames(res) <- c(paste("initialState.",stateInfo$genes,sep=""),
         paste("nextState.",stateInfo$genes,sep=""),
         "attractorAssignment","transitionsToAttractor")
  class(res) <- c("TransitionTable","data.frame")
  return(res)
}
