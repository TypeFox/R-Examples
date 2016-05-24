# Perform a Markov chain simulation with <numIterations> iterations/matrix multiplications on <network>.
# If <startStates> is supplied, probabilities for all other start states will be set to 0.
# All probabilities below <cutoff> are regarded as zero. 
markovSimulation <- function(network, numIterations=1000, startStates=list(), cutoff=0.001, returnTable=TRUE)
{
  stopifnot(inherits(network,"ProbabilisticBooleanNetwork") 
            | inherits(network,"BooleanNetwork"))

  if (sum(network$fixed == -1) > 32)
    stop("A Markov chain simulation with more than 32 non-fixed genes is not supported!")

  # the C code requires all interactions to be coded into one vector:

  if (inherits(network,"BooleanNetwork"))
  # deterministic network
  {
    # Assemble all input gene lists in one list <inputGenes>, and remember the split positions in <inputGenePositions>.
    inputGenes <- as.integer(unlist(lapply(network$interactions,function(interaction)interaction$input)))
    inputGenePositions <- as.integer(cumsum(c(0,sapply(network$interactions,
                                     function(interaction)length(interaction$input)))))

    # Do the same for the transition functions.
    transitionFunctions <- as.integer(unlist(lapply(network$interactions,function(interaction)interaction$func)))
    transitionFunctionPositions <- as.integer(cumsum(c(0,sapply(network$interactions,
                                              function(interaction)length(interaction$func)))))
    probabilities <- as.double(rep(1.0,length(network$genes)))
    functionAssignment <- as.integer(0:(length(network$genes)-1))
  }
  else
  # probabilistic network
  {
    wrongProb <- sapply(network$interactions,function(interaction)
                                    abs(1.0-sum(sapply(interaction,function(func)func$probability))) > 0.0001)
    if (any(wrongProb))
      stop(paste("The probabilities of gene(s) ",paste(network$genes[wrongProb],collapse=", ")," do not sum up to 1!",sep=""))

    # Assemble all input gene lists in one list <inputGenes>, and remember the split positions in <inputGenePositions>.
    inputGenes <- as.integer(unlist(lapply(network$interactions,function(interaction)lapply(interaction,function(singleFunc)singleFunc$input))))
    inputGenePositions <- as.integer(cumsum(c(0,unlist(lapply(network$interactions,
             function(interaction)lapply(interaction,function(singleFunc)length(singleFunc$input)))))))

    # Do the same for the transition functions.  
    transitionFunctions <- as.integer(unlist(lapply(network$interactions,function(interaction)lapply(interaction,function(singleFunc)singleFunc$func))))
  
    transitionFunctionPositions <- as.integer(cumsum(c(0,unlist(lapply(network$interactions,
             function(interaction)lapply(interaction,function(singleFunc)length(singleFunc$func)))))))

    probabilities <- as.double(unlist(lapply(network$interactions,function(interaction)lapply(interaction,function(singleFunc)singleFunc$probability))))

    functionAssignment <- as.integer(unlist(mapply(function(index,interaction)rep(index,length(interaction)),0:(length(network$interactions )- 1),
                              network$interactions)))
  }

  # check whether the start states are allowed
  # by comparing the values of the fixed genes
  fixedGenes <- which(network$fixed != -1)

  if (length(startStates) > 0)
  {
    statesValid <- sapply(startStates,function(state)
                      {
                        isTRUE(all(state[fixedGenes] == network$fixed[fixedGenes]))
                      })
    if (!isTRUE(all(statesValid)))
      warning("Some of the supplied start states did not match the restrictions of the fixed genes and were removed!")

    startStates <- startStates[statesValid]
  }

  convertedStartStates <- NULL

  if (length(startStates) > 0)
    convertedStartStates <- sapply(startStates,function(x)bin2dec(x,length(network$genes)))

  on.exit(.C("freeAllMemory", PACKAGE = "BoolNet"))
  # call C code
  res <- .Call("markovSimulation_R",
      inputGenes,inputGenePositions,
      transitionFunctions,transitionFunctionPositions,
      as.integer(network$fixed),
      functionAssignment,
      probabilities,
      as.integer(numIterations),
      convertedStartStates,
      as.double(cutoff),
      as.integer(returnTable),
      PACKAGE="BoolNet")

  if (length(network$genes) %% 32 == 0)
    numElementsPerEntry <- as.integer(length(network$genes) / 32)
  else
    numElementsPerEntry <- as.integer(length(network$genes) / 32  + 1)

  # build result matrix
  reachedStates <- data.frame(t(sapply(res$states,function(state)dec2bin(state,length(network$genes)))),res$probabilities)
  colnames(reachedStates) <- c(network$genes,"Probability") 
  
  if (returnTable)
  {
    initialStates <- matrix(res$initialStates,nrow=numElementsPerEntry)
    nextStates <- matrix(res$nextStates,nrow=numElementsPerEntry)
    res <- list(reachedStates=reachedStates,
                table=list(initialStates=initialStates,nextStates=nextStates,probabilities=res$transitionProbabilities),
                genes=network$genes)
  }
  else
  {
    res <- list(reachedStates=reachedStates,
                table=NULL,
                genes=network$genes)
  }
              
  class(res) <- "MarkovSimulation"
  return(res)
}
