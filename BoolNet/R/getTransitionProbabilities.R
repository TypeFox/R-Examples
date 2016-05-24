# Builds a matrix of transitions and their probabilities from <markovSimulation>
getTransitionProbabilities <- function(markovSimulation)
{
  stopifnot(inherits(markovSimulation,"MarkovSimulation"))
  if (is.null(markovSimulation$table))
    stop(paste("The supplied simulation result does not contain transition information.",
               "Please re-run markovSimulation() with returnTable=TRUE!"))
  
  initialStates <- t(apply(markovSimulation$table$initialStates,2,dec2bin,length(markovSimulation$genes)))

  nextStates <- t(apply(markovSimulation$table$nextStates,2,dec2bin,length(markovSimulation$genes)))

  idx <- order(apply(initialStates,1,function(x)paste(x,collapse="")))
  
  res <- data.frame(initialStates,nextStates,markovSimulation$table$probabilities)
  
  res <- res[idx,]
  
  rownames(res) <- NULL  
  colnames(res) <- c(paste("initialState.",markovSimulation$genes,sep=""),
                     paste("nextState.",markovSimulation$genes,sep=""),
                     "probability")
  return(res)
}
