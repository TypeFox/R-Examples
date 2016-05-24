simulateSymbolicModel <- function(network, 
                          method=c("exhaustive","random","chosen","sat.exhaustive","sat.restricted"),
						  startStates=NULL, 
						  returnSequences=(!(match.arg(method) %in% c("sat.exhaustive", "sat.restricted"))),
						  returnGraph=(!(match.arg(method) %in% c("sat.exhaustive", "sat.restricted"))),
						  returnAttractors=TRUE,
						  maxTransitions=Inf,
						  maxAttractorLength=Inf,
						  canonical=TRUE)
{
  stopifnot(inherits(network,"SymbolicBooleanNetwork"))
  if (is.null(network$internalStructs) || checkNullPointer(network$internalStructs))
  # refresh internal tree representations if necessary
    network$internalStructs = .Call("constructNetworkTrees_R",network);
  
  if (!is.null(maxAttractorLength) && is.infinite(maxAttractorLength))
    maxAttractorLength <- NULL
  
  if (length(method) > 1)
  # if no method is supplied, infer method from the type of <startStates>
  {
    if (!is.null(maxAttractorLength))
      method <- "sat.restricted"
    else
    if (length(startStates) == 0)
    {
      method <- "exhaustive"
    }
    else
    if (is.numeric(startStates))
    {
      if (length(startStates) > 1)
        stop("Please supply either the number of start states or a list of start states in startStates!")
      else
        method <- "random"
    }
    else
    if (is.list(startStates))
      method <- "chosen"
    else
      stop("Please supply either the number of start states or a list of start states in startStates!")
  }
  
  method <- match.arg(method, c("exhaustive","random","chosen","sat.exhaustive","sat.restricted"))
  
  if (method == "random")
  {
    if (startStates > 2 ^ (sum(network$timeDelays[network$fixed == -1])))
          # more start states than in the full network
     {
              method <- "exhaustive"
              warning(paste("The number of random states is set larger than the total",
                            "number of states. Performing an exhaustive search!"))
     }                            
  }
    
  if (method %in% c("sat.exhaustive", "sat.restricted"))
  {
    if (!is.null(maxAttractorLength))
      maxAttractorLength <- as.integer(maxAttractorLength)
    else
	if (method == "sat.restricted")
		stop("maxAttractorLength must be set for method=\"sat.restricted\"!")
  
    if (returnSequences)
    {
      warning("Sequences cannot be returned for method=\"sat.exhaustive\" and method=\"sat.restricted\"!")
      returnSequences <- FALSE
    }
    if (returnGraph)
    {
      warning("Graph cannot be returned for method=\"sat.exhaustive\" and method=\"sat.restricted\"!")
      returnGraph <- FALSE
    }
  }
  else
  if (method == "exhaustive")
  {
    startStates <- NULL
    convertedStates <- NULL
  }
  else
  if (method == "chosen")
  {
    convertedStates <- 
      lapply(startStates, function(state)
      {
        if (!is.null(dim(state)))
        {
          if (ncol(state) != length(network$genes))
            stop(paste("\"startStates\" must be either a list of vectors with one value for each gene,",
                       "or a list of matrices with the genes in the columns and multiple predecessor states in the rows!"))
          return(as.integer(as.matrix(state[nrow(state):1,])))
        }
        else
        if (!is.numeric(state) || length(state) != length(network$genes))
          stop(paste("\"startStates\" must be either a list of vectors with one value for each gene,",
                     "or a list of matrices with the genes in the columns and multiple predecessor states in the rows!"))
         return(as.integer(state))
      })
  }
  else
    convertedStates <- as.integer(startStates)
  
  if (maxTransitions == 0)
    warning("\"maxTransitions\" is set to 0, which disables the transition limit!")
  else
  if (is.infinite(maxTransitions))
    maxTransitions <- 0
    
  on.exit(.C("freeAllMemory", PACKAGE = "BoolNet"))
  
  if (method %in% c("sat.exhaustive","sat.restricted"))  
    res <- .Call("symbolicSATSearch_R", network$internalStructs, maxAttractorLength, method == "sat.restricted")
  else
    res <- .Call("simulateStates_R",
                 network$internalStructs, 
                 convertedStates, 
                 #as.integer(length(startStates)), 
                 as.integer(maxTransitions),
                 as.logical(returnSequences),
                 as.logical(returnGraph),
                 as.logical(returnAttractors))
 
  ret <- list()
  if (returnSequences)
  { 
    ret[["sequences"]] <- list()
  }
  
  if (returnGraph)
  {
    initialStates <- matrix(res[[2]][[1]], ncol=length(network$genes), byrow=TRUE)
    colnames(initialStates) <- paste("initialState.",network$genes, sep="")
    
    nextStates <- matrix(res[[2]][[2]], ncol=length(network$genes), byrow=TRUE)
    colnames(nextStates) <- paste("nextState.",network$genes, sep="")
    
    attractorAssignment <- data.frame(res[[2]][[3]] + 1)
    attractorAssignment[attractorAssignment == 0] <- NA
    colnames(attractorAssignment) <- "attractorAssignment"
    graph <- data.frame(unique(cbind(initialStates, nextStates, attractorAssignment)))
    class(graph) <- c("TransitionTable","data.frame")
    ret[["graph"]] <- graph
  }
  
  if (returnAttractors)
  {
    attractors <- lapply(res[[3]], function(x)
    {
      att <- matrix(x, ncol=length(network$genes), byrow=TRUE)
      colnames(att) <- network$genes
      
      if (canonical)
      {
        smallestIndex <- -1
        smallestVal  <- rep(Inf, ncol(att))
        for (i in seq_len(nrow(att)))
        # iterate over elements of encoded state
        {
		  equal <- TRUE
          for (j in seq(ncol(att),by=-1,length.out=ncol(att)))
          {
            if (att[i,j] < smallestVal[j])
            # determine new minimum
            {
			  equal <- FALSE				
              smallestVal <- att[i,]
              smallestIndex <- i
              break
            }
            else
			if (att[i,j] > smallestVal[j])
			{
				equal <- FALSE        
				break
			}
          }
		  if (equal && i != smallestIndex)
		  {
			  if (i - smallestIndex < (smallestIndex + nrow(att) - i) %% nrow(att))
			  {
				  smallestVal <- att[i,]
				  smallestIndex <- i
			  }
		  }
        }
      }
      if (smallestIndex != 1)
        # rearrange matrix
        att <- rbind(att[smallestIndex:nrow(att),,drop=FALSE],
                     att[seq_len(smallestIndex-1),,drop=FALSE])
      
      as.data.frame(att)
    })
    ret[["attractors"]] <- attractors
  }
  
  if (returnSequences)
  {               
    maxDelay <- max(network$timeDelays)
    
    if (returnAttractors)
    {
      ret[["attractorAssignment"]] <- res[[4]] + 1
      ret[["attractorAssignment"]][ret[["attractorAssignment"]] == 0] <- NA
    }
    
    sequences <- mapply(function(seq, att)
    {
      seq <- matrix(seq, ncol=length(network$genes), byrow=TRUE)
      colnames(seq) <- network$genes
      rownames(seq) <- paste("t =", (1:nrow(seq)) - maxDelay)
      seq <- as.data.frame(seq)
      
      # add attractor information to sequence
      if (!is.na(att))
        attributes(seq)$attractor <-
          (nrow(seq) - nrow(ret[["attractors"]][[att]]) + 1):nrow(seq)
      seq
    }, res[[1]], ret[["attractorAssignment"]], SIMPLIFY=FALSE)
    ret[["sequences"]] <- sequences
  }
  class(ret) <- "SymbolicSimulation"  
  return(ret)
}
