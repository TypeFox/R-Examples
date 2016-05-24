
perturbTrajectories <- function(network, 
                                measure = c("hamming", "sensitivity", "attractor"),
                                numSamples = 1000,
                                flipBits = 1,
                                updateType = c("synchronous","asynchronous","probabilistic"),
                                gene,
                                ...)
{
  symbolic <- inherits(network,"SymbolicBooleanNetwork")
  probabilistic <- inherits(network,"ProbabilisticBooleanNetwork") 
  stopifnot(inherits(network,"BooleanNetwork") || symbolic || probabilistic)
  
  measure <- match.arg(measure, c("hamming", "sensitivity", "attractor"))

  if (probabilistic && measure != "hamming")
    stop("Probabilistic networks are only supported for measure=\"hamming\"!")

  fixedIdx <- which(network$fixed != -1)
  nonFixedIdx <- which(network$fixed == -1)  

  if (symbolic)
  {
    maxDelay <- max(network$timeDelays)
    flipIndices <- unlist(mapply(function(x, i)
                          {
                           if (network$fixed[i] == -1)
                             seq_len(x) + (i - 1) * maxDelay
                           else
                            c()  
                          }, network$timeDelays, seq_along(network$genes)))
    startStates <- lapply(seq_len(numSamples), function(i)
    {
      m <- matrix(nrow=maxDelay, ncol=length(network$genes), 
             round(runif(n=maxDelay*length(network$genes))))
      m[,fixedIdx] <- network$fixed[fixedIdx]
      m
    })
  }
  else
  {
    flipIndices <- seq_along(network$genes)[nonFixedIdx]
    startStates <- lapply(seq_len(numSamples), function(i)
    {
      s <- round(runif(n=length(network$genes)))
      s[fixedIdx] <- network$fixed[fixedIdx]
      s
    })
  }
  
  if (length(flipIndices) == 0)
    stop("All genes in this network are fixed!")
  
  perturbedStates <- lapply(startStates, function(state)
  {
    flip <- sample(flipIndices, size=flipBits, replace=FALSE)
    state[flip] <- 1 - state[flip]
    return(state)
  })
  

  if (measure == "hamming")
  {
    dists <- mapply(function(state1, state2)
                     {
                       sum(stateTransition(network=network,
                                           state=state1,
                                           type=updateType,
                                           ...) !=
                           stateTransition(network=network,
                                           state=state2,
                                           type=updateType,
                                           ...))                                          
                     }, startStates, perturbedStates)/length(network$genes)
    return(list(stat=dists, value=mean(dists)))                    
  }
  else
  if (measure == "sensitivity")
  {
    if (missing(gene))
      stop("Parameter \"gene\" must be set for sensitivity analysis!")
    
    if (is.character(gene))
      geneIdx <- which(network$genes == gene)
    else
      geneIdx <- gene
    
    if (length(geneIdx) == 0)
      stop(paste("Unknown gene \"", gene, "\"!", sep=""))
    
    sensitivities <- unname(mapply(function(state1, state2)
                     {
                       stateTransition(network=network,
                                           state=state1,
                                           chosenGene = geneIdx,
                                           ...)[geneIdx] !=
                       stateTransition(network=network,
                                           state=state2,
                                           chosenGene = geneIdx,
                                           ...)[geneIdx]                                          
                     }, startStates, perturbedStates))
    return(list(stat=sensitivities, value=mean(sensitivities)))
  }
  else
  {
    if (symbolic)
    {
      att1 <-  simulateSymbolicModel(network, 
                                     startStates=startStates, 
                                     returnGraph=FALSE,
                                     canonical=TRUE,                                      
                                     ...)
      att2 <- simulateSymbolicModel(network, 
                                    startStates=perturbedStates, 
                                    returnGraph=FALSE, 
                                    canonical=TRUE,                                    
                                    ...)
      attractorAssignment1 <- att1$attractorAssignment
      attractorAssignment2 <- att2$attractorAssignment    
    }
    else
    {
      att1 <- getAttractors(network, 
                            startStates=startStates, 
                            returnTable=TRUE,
                            canonical=TRUE,                                         
                            ...)
      att2 <- getAttractors(network, 
                            startStates=perturbedStates, 
                            returnTable=TRUE,
                            canonical=TRUE,                            
                             ...)
                             
      tt1 <- getTransitionTable(att1)
      tt2 <- getTransitionTable(att2)
      attractorAssignment1 <- tt1$attractorAssignment
      attractorAssignment2 <- tt2$attractorAssignment
      names(attractorAssignment1) <- apply(tt1[,seq_along(network$genes)],1,paste,collapse="")
      names(attractorAssignment2) <- apply(tt2[,seq_along(network$genes)],1,paste,collapse="")
      
      attractorAssignment1 <- unname(attractorAssignment1[sapply(startStates,paste,collapse="")])
      attractorAssignment2 <- unname(attractorAssignment2[sapply(perturbedStates,paste,collapse="")])
      
    }
    
    attractorsEqual <- matrix(sapply(seq_along(att1$attractors), function(i1)
                       {
                         seq1 <- getAttractorSequence(att1, i1)
                         sapply(seq_along(att2$attractors), function(i2)
                         {
                            seq2 <- getAttractorSequence(att2, i2)
                            return(identical(seq1,seq2))
                         })
                       }), nrow=length(att1$attractors))
                           
    stat <- apply(cbind(attractorAssignment1, attractorAssignment2), 1, function(idx)
            {
              attractorsEqual[idx[1],idx[2]]
            })
    return(list(stat=stat, value=mean(stat)))
  }          
}
