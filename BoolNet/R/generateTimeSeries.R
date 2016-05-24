
# Generate time series from a Boolean network <network>.
# <numSeries> is the number of start states for which time series are generated.
# <numMeasurements> is the number of time points for each time series.
# <type> is the type of transitions used (synchronous or asynchronous)
# If <type> is "asynchronous", <geneProbabilities> describes the 
# transition probabilities for the genes.
# If <noiseLevel> is not 0, Gaussian noise is added to the result.
generateTimeSeries <- function(network, numSeries, numMeasurements, 
                            type = c("synchronous","asynchronous","probabilistic"),
                            geneProbabilities, 
                            perturbations=0,
                            noiseLevel = 0.0)
{
  symbolic <- inherits(network, "SymbolicBooleanNetwork")

  if (missing(geneProbabilities))
    geneProbabilities <- NULL
  
  if (perturbations > 0)
  {
    perturbationMatrix <- sapply(1:numSeries, function(x)
                                 {
                                   p <- rep(NA, length(network$genes))
                                   p[sample(1:length(network$genes),size=perturbations)] <- 
                                      sample(c(0,1),size=perturbations,replace=TRUE)
                                   
                                   return(p)
                                 })
    rownames(perturbationMatrix) <- network$genes                                 
  }
    
  ts <- lapply(seq_len(numSeries), function(i)
  {
    if (symbolic)
    {
      if (perturbations > 0)
      {
        fixedIdx <- which(perturbationMatrix[,i] != -1)
        network <- fixGenes(network, fixedIdx, perturbationMatrix[fixedIdx,i])
      }
      
      res <- t(simulateSymbolicModel(network, startStates=1, 
                                     returnAttractors=FALSE, 
                                     returnGraph=TRUE, 
                                     returnSequences=TRUE)$sequences[[1]])
    }
    else
    {
      startState <- round(runif(length(network$genes)))
      
      if (perturbations > 0)
      {
        fixedIdx <- which(perturbationMatrix[,i] != -1)
        network <- fixGenes(network, fixedIdx, perturbationMatrix[fixedIdx,i])
        startState[fixedIdx] <- perturbationMatrix[fixedIdx,i]
      }
      
      res <- startState
      for (j in 2:numMeasurements)
      {
        startState <- stateTransition(network, startState, 
                                      type=type, geneProbabilities=geneProbabilities)
        res <- cbind(res,startState)
      }
    }    
    colnames(res) <- NULL

    if (noiseLevel != 0)
    {
      res <- res + matrix(rnorm(mean=0, sd=noiseLevel, n = length(res)), nrow=nrow(res))
    }
    
    rownames(res) <- network$genes
    colnames(res) <- seq_len(ncol(res))
    return(res)
  })
  
  if (perturbations > 0)
    ts$perturbations <- perturbationMatrix
  
  return(ts)
}
