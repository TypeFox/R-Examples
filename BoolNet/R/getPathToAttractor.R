getPathToAttractor <- function(network, state, includeAttractorStates=c("all","first","none"))
{
  stopifnot(inherits(network,"BooleanNetwork") || inherits(network,"SymbolicBooleanNetwork") 
            || inherits(network, "AttractorInfo"))
  
  includeAttractorStates <- match.arg(includeAttractorStates, c("all","first","none"))
  
  if (inherits(network,"SymbolicBooleanNetwork"))
  {
    sim <- simulateSymbolicModel(network, startStates=list(state))
    res <- data.frame(sim$sequences[[1]])
    
    attractorStates <- nrow(res) - seq(length.out=nrow(sim$attractors[[1]]),by=-1,to=1) + 1
    
    if (includeAttractorStates == "first")
    {
      res <- res[seq_len(attractorStates[1]),,drop=FALSE]
      attractorStates <- attractorStates[1]
    }
    else
    if (includeAttractorStates == "none")
    {
      res <- res[-attractorStates,,drop=FALSE]
      attractorStates <- NULL
    }
    attributes(res)$attractor <- attractorStates
  }
  else
  {  
    if (inherits(network,"BooleanNetwork"))
    {
      table <- getTransitionTable(getAttractors(network, startStates=list(state)))
    }
    else
    {
      if (is.null(network$stateInfo$table))
        stop(paste("This AttractorInfo structure does not contain transition table information.",
             "Please re-run getAttractors() with a synchronous search and returnTable=TRUE!"))
      table <- getTransitionTable(network)
    }    
    
    numGenes <- (ncol(table) - 2) / 2
    initialStates <- apply(table[,seq_len(numGenes),drop=FALSE],1,function(x)paste(x,collapse=""))
    
    currentState <- state
    
    res <- data.frame(matrix(state,nrow=1))
    
    attractorStart <- NA
    stateCount <- 1
    repeat
    {
      currentStateIdx <- which(initialStates == paste(currentState,collapse=""))

      if (length(currentStateIdx) == 0)
        stop(paste("Could not find state",paste(currentState,collapse=""),"in the transition table!"))

      if (table[currentStateIdx,"transitionsToAttractor"] == 0 && is.na(attractorStart))
        attractorStart <- stateCount

      # stop depending on "includeAttractorStates" option
      if ((includeAttractorStates == "all" && stateCount == nrow(table))
          || (includeAttractorStates == "first" && table[currentStateIdx,"transitionsToAttractor"] == 0)
          || (includeAttractorStates == "none" && table[currentStateIdx,"transitionsToAttractor"] <= 1))
        break
        
      currentState <- as.integer(table[currentStateIdx,(numGenes+1):(2*numGenes)])
      res <- rbind(res,currentState)
      
      stateCount <- stateCount + 1
    }
    
    if (!is.na(attractorStart))
      attractorIdx <- seq(attractorStart, nrow(table), by=1)
    else
      attractorIdx <- NULL
    
    # special case: start state is attractor state and we do not want to include attractor states
    # => return empty data frame
    if (includeAttractorStates == "none" && table[currentStateIdx,"transitionsToAttractor"] == 0)
    {
      res <- data.frame(matrix(nrow=0,ncol=numGenes))
    }
    else
      if (length(attractorIdx) > 0)
        attributes(res)$attractor <- attractorIdx

    colnames(res) <- sapply(colnames(table)[seq_len(numGenes)],function(n)strsplit(n,".",fixed=TRUE)[[1]][2])
    
    rownames(res) <- NULL
  }
  return(res)
}
