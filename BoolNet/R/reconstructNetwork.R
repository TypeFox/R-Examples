# Reconstruct a Boolean network from a transition table or a list of time series in <measurements>.
# If <method> is "bestfit", Lähdesmäki's best-fit extension algorithm is called.
# If <method> is "reveal", Liang's REVEAL algorithm is called.
# <maxK> specifies the maximum number of input genes of a function.
# If <readableFunctions> is true, DNF representations of the functions are generated.
reconstructNetwork <- function(measurements,method=c("bestfit","reveal"),maxK=5,
                               requiredDependencies = NULL, excludedDependencies = NULL, perturbations=NULL,
                               readableFunctions=FALSE, allSolutions=FALSE, returnPBN=FALSE)
{
  if (maxK < 0)
    stop("maxK must be >= 0!")

  # determine method to use
  meth <- switch(match.arg(method,c("bestfit","reveal")),
      bestfit=0,
      reveal=1,
      stop("'method' must be one of \"bestfit\",\"reveal\""))
      
  perturbationMatrix <- NULL     
      
  if (inherits(measurements,"TransitionTable"))
  # preprocess transition table and call algorithm
  {
    numGenes <- (ncol(measurements) - 2) / 2
    
    if (numGenes < maxK)
    {
      maxK <- numGenes
      warning(paste("maxK was chosen greater than the total number of input genes and reset to ",numGenes,"!",sep=""))
    }
    
    if (!is.null(perturbations))
    { 
      if (!all(perturbations %in% c(0,1,-1,NA)))
        stop("The perturbation matrix may only contain the values 0,1,-1 and NA!")
      
      perturbations[is.na(perturbations)] <- -1 
        
      if (is.null(dim(perturbations)))
        perturbations <- matrix(perturbations,ncol=1)
        
      if (ncol(perturbations) != 1 && ncol(perturbations) != nrow(measurements))
        stop(paste("There must be either one global vector of perturbations, or a matrix",
                   "containing exact one column for each time series!"))
     
      if (nrow(perturbations) != numGenes)
        stop("The perturbation matrix must have exactly the same number of rows as the measurements!")
        
      if (ncol(perturbations) == 1)
      {
        perturbationMatrix <- matrix(-1,nrow=numGenes,ncol=nrow(measurements))
        for (j in seq_len(nrow(measurements)))
          perturbationMatrix[,j] <- perturbations[,1]
      }
      else
        perturbationMatrix <- as.matrix(perturbations)
    }   
    
    inputStates <- as.integer(t(as.matrix(measurements[,seq_len(numGenes)])))
    outputStates <- as.integer(t(as.matrix(measurements[,(numGenes+1):(2*numGenes)])))
    numStates <- nrow(measurements)
    genenames <- sapply(colnames(measurements)[seq_len(numGenes)],function(x)strsplit(x,".",fixed=TRUE)[[1]][2])    
  }
  else
  # the measurements are time series
  {
    if (!is.null(dim(measurements)))
    # only one time series => create list
      measurements <- list(measurements)
      
    numGenes <- nrow(measurements[[1]])
    
    if (numGenes < maxK)
    {
      maxK <- numGenes
      warning(paste("maxK was chosen greater than the total number of input genes and reset to ",numGenes,"!",sep=""))
    }

    if (is.null(perturbations) && !is.null(measurements$perturbations))
    {
      perturbations <- measurements$perturbations
      measurements$perturbations <- NULL
    }
    
    if (!is.null(perturbations))
    { 
      if (!all(perturbations %in% c(0,1,-1,NA)))
        stop("The perturbation matrix may only contain the values 0,1,-1 and NA!")
      
      perturbations[is.na(perturbations)] <- -1    
      
      if (is.null(dim(perturbations)))
        perturbations <- matrix(perturbations,ncol=1)

      if (ncol(perturbations) != 1 && ncol(perturbations) != length(measurements))
        stop(paste("There must be either one global vector of perturbations, or a matrix",
                   "containing exact one column for each time series!"))                   
     
      if (nrow(perturbations) != numGenes)
        stop("The perturbation matrix must have exactly the same number of rows as the measurements!")
    }
    
    perturbationMatrix <- c()
    
    genenames <- rownames(measurements[[1]])
    if (is.null(genenames))
      genenames <- paste("Gene",seq_len(numGenes))
    inputStates <- c()
    outputStates <- c()
    
    for (i in seq_along(measurements))
    # iterate over all time series and build state vectors
    {
      measurement <- measurements[[i]]
      if (numGenes != nrow(measurement))
        stop("All measurement matrices must contain the same genes!")
      inputStates <- c(inputStates,as.integer(as.matrix(measurement[,1:(ncol(measurement)-1)])))
      outputStates <- c(outputStates,as.integer(as.matrix(measurement[,2:ncol(measurement)])))

      if (!is.null(perturbations))
      {
        for (j in seq_len(ncol(measurement) - 1))
        {
          if (ncol(perturbations) == 1)
            perturbationMatrix <- cbind(perturbationMatrix, perturbations[,1])
          else
            perturbationMatrix <- cbind(perturbationMatrix, perturbations[,i])
        }
      }
    }
    numStates <- as.integer(length(inputStates) / numGenes)
  }
  
  if (is.null(requiredDependencies))
    requiredDepMatrix <- NULL
  else
  {
    if (is.null(names(requiredDependencies)))
      names(requiredDependencies) <- genenames
    
    if (length(union(names(requiredDependencies),genenames)) != length(genenames))
      stop("The required dependencies must consist of gene names that are comprised in the measurements!")

    requiredDepMatrix <- matrix(0,nrow=numGenes,ncol=numGenes)
    colnames(requiredDepMatrix) <- genenames
    rownames(requiredDepMatrix) <- genenames
    
    maxRequired <- 0      
    for (i in seq_along(requiredDependencies))
    {
      maxRequired <- max(maxRequired, length(requiredDependencies[i]))
      for (el in requiredDependencies[i])
        requiredDepMatrix[el,names(requiredDependencies)[i]] <- 1
    }
    
    if (maxRequired > maxK)
    {
      warning(paste("The number of required dependencies is greater than maxK! Setting maxK to ",                  maxRequired, "!", sep=""))
      maxK <- maxRequired
    } 
    
    requiredDepMatrix <- as.integer(requiredDepMatrix) 
  }
  
  if (is.null(excludedDependencies))
    excludedDepMatrix <- NULL
  else
  {
    if (is.null(names(excludedDependencies)))
      names(excludedDependencies) <- genenames
    
    if (length(union(names(excludedDependencies),genenames)) != length(genenames))
      stop("The excluded dependencies must consist of gene names that are comprised in the measurements!")

    excludedDepMatrix <- matrix(0,nrow=numGenes,ncol=numGenes)
    colnames(excludedDepMatrix) <- genenames
    rownames(excludedDepMatrix) <- genenames
          
    for (i in seq_along(excludedDependencies))
    {
      gene <- names(excludedDependencies)[i]
      if (!is.null(requiredDependencies) && 
          !is.null(requiredDependencies[[gene]]))
      {
        conflicts <- intersect(requiredDependencies[[gene]],
                               excludedDependencies[[gene]])
        
        if (length(conflicts) > 0)
          stop(paste("For gene ",gene,
                      ", potential inputs were specified both as required and excluded dependencies: ",
                      paste(conflicts, collapse=", "), sep=""))
      }
      for (el in excludedDependencies[i])
        excludedDepMatrix[el,names(excludedDependencies)[i]] <- 1
    }
    excludedDepMatrix <- as.integer(excludedDepMatrix)
  }  
  
  if (!is.null(perturbations))
    perturbationMatrix <- as.integer(perturbationMatrix)
  
  on.exit(.C("freeAllMemory", PACKAGE = "BoolNet"))
  # call C code
  res <- .Call("reconstructNetwork_R",
      inputStates,
      outputStates,
      perturbationMatrix,
      as.integer(numStates),
      requiredDepMatrix,
      excludedDepMatrix,
      as.integer(maxK),
      as.integer(meth),
      as.integer(allSolutions),
      as.integer(returnPBN))
  
  
  if (any(sapply(res,function(interaction)length(interaction)==0)))
  # some function lists are empty
    warning("Some functions could not be inferred. Possibly the input data is noisy or maxK was chosen too small!")
    
  # prepare result object
  res <- list(genes=genenames,
              interactions=lapply(res,function(gene)
                    lapply(gene,function(interaction)
                    {
                      interaction$expression <- 
                                getInteractionString(readableFunctions,
                                           interaction$func,
                                           genenames[interaction$input])
                      if (returnPBN)
                        interaction$probability <- 1.0/length(gene)
                        
                      interaction
                    })),
              fixed=sapply(res,function(gene)
                  {
                    if (length(gene) == 0)
                      -1
                    else
                    if (gene[[1]]$input[1] == 0)
                      gene[[1]]$func[1]
                    else
                      -1
                  }))
  
  names(res$interactions) <- res$genes
  names(res$fixed) <- res$genes
  
  if (returnPBN)
    class(res) <- "ProbabilisticBooleanNetwork"
  else
    class(res) <- "BooleanNetworkCollection"
  
  
  if (allSolutions)
  # simplify functions and remove duplicates
  {
    res <- simplifyNetwork(res)
    res$interactions <- lapply(res$interactions,function(interaction)
                               {
                                 duplicates <- duplicated(sapply(interaction,function(func)func$expression))
                                 return(interaction[!duplicates])
                               })
  }
  
  return(res)
}
