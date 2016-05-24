# Summarize information on a state <state>, specified as a vector of
# Boolean values. Information is taken from the state tables in <attractorInfo>.
# The function outputs the corresponding attractor, the distance to the basin, 
# and the next state.
getStateSummary <- function(attractorInfo,state)
{
  stopifnot(inherits(attractorInfo,"AttractorInfo") | inherits(attractorInfo,"SymbolicSimulation"))
  
  if (inherits(attractorInfo,"SymbolicSimulation"))
  {
    if (is.null(attractorInfo$graph))
      stop(paste("This SymbolicSimulation structure does not contain transition table information.",
           "Please re-run simulateSymbolicModel() with returnGraph=TRUE!"))
    
    geneCols <- setdiff(colnames(attractorInfo$graph),c("attractorAssignment","transitionsToAttractor"))
    numGenes <- (length(geneCols)) / 2
    
    stateIndices <- apply(attractorInfo$graph[,seq_len(numGenes),drop=FALSE], 1, function(x)
                    {
                      all(as.integer(x) == as.integer(state))
                    })
           
    return(attractorInfo$graph[stateIndices,,drop=FALSE])
  }
  else
  {
    if (length(state) != length(attractorInfo$stateInfo$genes))
      stop("State must have one element for each gene!")
    
    if (is.null(attractorInfo$stateInfo$table))
      stop(paste("This AttractorInfo structure does not contain transition table information.",
             "Please re-run getAttractors() with a synchronous search and returnTable=TRUE!"))

    if (!is.null(attractorInfo$stateInfo$initialStates))
    {
      stateNo <- bin2dec(state,length(state))
      
      # search element in initial state index
      stateNo <- which(apply(attractorInfo$stateInfo$initialStates,2,function(x)
          {
            isTRUE(all.equal(x,stateNo))
          }))
    }
    else
    {
      # find out the decimal number of the state which is used in the
      # transition table
      shortenedState <- state[attractorInfo$stateInfo$fixedGenes == -1]
      
      # simply use <stateNo>-th elements if initial states are not provided    
      stateNo <- bin2dec(shortenedState,length(shortenedState)) + 1
    }
    
    # coerce summary into a one-row dataframe, i.e. create a
    # TransitionTable object with one element
    res <- as.data.frame(as.list(unlist(list(state,
          dec2bin(attractorInfo$stateInfo$table[,stateNo],
            length(attractorInfo$stateInfo$genes)),
          attractorInfo$stateInfo$attractorAssignment[stateNo],
          attractorInfo$stateInfo$stepsToAttractor[stateNo]))))
          
    colnames(res) <- c(paste("initialState.",attractorInfo$stateInfo$genes,sep=""),
           paste("nextState.",attractorInfo$stateInfo$genes,sep=""),
           "attractorAssignment","transitionsToAttractor")  
    class(res) <- c("TransitionTable","data.frame")
    
    return(res)
  }  
}
