# Identify attractors in a Boolean network.
# <network> is a BooleanNetwork/SymbolicBooleanNetwork structure specifying the network.
# <method> specifies what kind of search is conducted: "exhaustive" performs a full search over all states,
# "random" generates <startStates> random initial states, and "chosen" uses the states supplied in <startStates>. 
# "sat.exhaustive" and "sat.restricted" start a SAT-based attractor search.
# <genesON> and <genesOFF> are lists of genes to be set to ON/1 or OFF/0 respectively.
# If <canonical> is true, states in the attractors are reordered such that the "smallest" state is the first
# <randomChainLength> is the number of random transitions performed for the identification of an asynchronous attractor
# If <avoidSelfLoops> is true, loops to the same state are eliminated from asynchronous attractors.
# <geneProbabilities> optionally specifies the probabilities of choosing a gene for an asynchronous update.
# <maxAttractorLength> specifies the maximum attractor length for method="sat.restricted" and the initial search length for method="sat.exhaustive".
# if <returnTable> is true, the transition table is included in the result.
getAttractors <- function (network, type=c("synchronous","asynchronous"), 
         method=c("exhaustive","sat.exhaustive","sat.restricted","random","chosen"), startStates=list(),
         genesON = c(), genesOFF = c(), canonical=TRUE,
         randomChainLength = 10000, avoidSelfLoops = TRUE, 
         geneProbabilities = NULL, 
		 maxAttractorLength=Inf, 
         returnTable=TRUE) 
{
  stopifnot(inherits(network,"BooleanNetwork") || inherits(network,"SymbolicBooleanNetwork"))
  
  symbolic <- inherits(network,"SymbolicBooleanNetwork")
  
  nonFixedPositions <- which(network$fixed == -1)
  
  type <- match.arg(type, c("synchronous","asynchronous"))
  
  if (type == "asynchronous")
  {
    if (symbolic)
      stop("Only synchronous updates are allowed for symbolic networks!")
    
    if (length(method) == 1 & match.arg(method) == "exhaustive")
      stop("Asynchronous attractor search cannot be performed in exhaustive search mode!")
    
    if (length(geneProbabilities) > 0 )
    {
      if (length(geneProbabilities) != length(network$genes))
        stop("Please supply exactly one probability for each gene!")
      if (abs(1.0 - sum(geneProbabilities)) > 0.0001)
        stop("The supplied gene probabilities do not sum up to 1!")
    }
    
  }
  
  if (!is.null(maxAttractorLength) && is.infinite(maxAttractorLength))
    maxAttractorLength <- NULL
  
  if (length(method) > 1)
  # if no method is supplied, infer method from the type of <startStates>
  {
    if (type == "asynchronous" & length(startStates) == 0)
    {
      startStates <- max(round(2 ^ sum(network$fixed == -1) / 20), 5)
      method = "random"
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
    if (is.list(startStates) & (length(startStates) > 0))
      method <- "chosen"
    else
    if (!is.null(maxAttractorLength))
      method <- "sat.restricted"
    else
      method <- "exhaustive"
  }
 
  # fix genes according to genesON and genesOFF
  if (length(genesON) > 0) 
  {
    network <- fixGenes(network,genesON,1)
  }
  if (length(genesOFF) > 0) 
  {
    network <- fixGenes(network,genesOFF,0)
  }
    
  if (symbolic)
  {
    return(simulateSymbolicModel(network, 
                                 method=method, 
                                 startStates=startStates,
                                 maxAttractorLength=maxAttractorLength,
                                 returnGraph=returnTable && !(match.arg(method) %in% c("sat.exhaustive", "sat.restricted")),
                                 returnSequences=FALSE,
                                 returnAttractors=TRUE,
                                 canonical=canonical))
  }
  else
  {
    method <- match.arg(method,c("exhaustive","sat.exhaustive","sat.restricted","random","chosen"))

	if (method == "sat.restricted" && is.null(maxAttractorLength))
		stop("maxAttractorLength must be set for method=\"sat.restricted\"!")
	
    if ((length(network$genes) > 29) && (method == "exhaustive") && (type == "synchronous"))
    {
      method <- "sat.exhaustive"
      warning("An exhaustive state space search is restricted to networks with at most 29 genes. Switching to the SAT-based exhaustive search, which supports more genes, but does not return a transition table!")
    }
    else
    if (method %in% c("sat.exhaustive", "sat.restricted") && type != "synchronous")
      stop("SAT-based search can only be used for synchronous networks!")
    
    startStates <- switch(method,
      exhaustive = list(),
      sat.exhaustive = list(),
	  sat.restricted = list(),
      random = {
          if (!is.numeric(startStates))
            stop("Please supply the number of random states in startStates!")

          if (startStates > (2 ^ length(nonFixedPositions)))
          # more start states than in the full network
          {
            if (type == "synchronous")
            {
              list()
              warning(paste("The number of random states is set larger than the total",
                            "number of states. Performing an exhaustive search!"))
            }
            else
            {
              warning(paste("The number of random states is set larger than the total ",
                          "number of states! The maximum number of different states is ",2 ^ length(nonFixedPositions)),"!",sep="")
              startStates = 2 ^ length(nonFixedPositions)            
            }     
          }
          # generate random matrix
          generateRandomStartStates(network, startStates)

         },
      chosen = {
          if (!is.list(startStates) | length(startStates) == 0)
            stop("No start states supplied!")
          if (!all(sapply(startStates,function(x)(length(x) == length(network$genes)))))
            stop(paste("Please provide binary vectors with",length(network$genes),
              "elements in startStates!"))
    
            fixedGenes <- which(network$fixed != -1)
            statesValid <- sapply(startStates,function(state)
                                  {
                                    isTRUE(all(state[fixedGenes] == network$fixed[fixedGenes]))
                                  })
        startStates <- startStates[statesValid]
            if (!any(statesValid))
              stop("None of the supplied start states matched the restrictions of the fixed genes!")
            if (!all(statesValid))
              warning("Some of the supplied start states did not match the restrictions of the fixed genes and were removed!")    
              
          startStates
         }
    )
    
    if (!is.null(maxAttractorLength))
    {
      if (!(method %in% c("sat.exhaustive", "sat.restricted")))
        stop("maxAttractorLength can only be used with method=\"sat.exhaustive\" or method=\"sat-res.ricted\"!")
      maxAttractorLength <- as.integer(maxAttractorLength)
    }
    specialInitialization <- NULL
    
    convertedStartStates <- NULL

    if (length(startStates) > 0)
      convertedStartStates <- sapply(startStates,function(x)bin2dec(x,length(network$genes)))

    # the C code requires all interactions to be coded into one vector:
    # Assemble all input gene lists in one list <inputGenes>, and remember the split positions in <inputGenePositions>.
    inputGenes <- as.integer(unlist(lapply(network$interactions,function(interaction)interaction$input)))
    inputGenePositions <- as.integer(cumsum(c(0,sapply(network$interactions,
             function(interaction)length(interaction$input)))))

    # Do the same for the transition functions.
    transitionFunctions <- as.integer(unlist(lapply(network$interactions,function(interaction)interaction$func)))
    transitionFunctionPositions <- as.integer(cumsum(c(0,sapply(network$interactions,
                function(interaction)length(interaction$func)))))

   searchType <- switch(type,
      synchronous = if (method == "sat.exhaustive") 2 else if (method == "sat.restricted") 3 else 0,
      asynchronous = 1)
    
    on.exit(.C("freeAllMemory", PACKAGE = "BoolNet"))
    # Call the C code
    result <- .Call("getAttractors_R",inputGenes,inputGenePositions,
          transitionFunctions,transitionFunctionPositions,
          as.integer(network$fixed),
          as.integer(convertedStartStates),
          as.integer(searchType),
          as.double(geneProbabilities),
          as.integer(randomChainLength),
          as.integer(avoidSelfLoops),
          as.integer(returnTable),
          maxAttractorLength,
          PACKAGE="BoolNet")
    
    if (is.null(result))
      stop("An error occurred in external C code!")
    
    if (length(result$attractors) == 0)
      stop("getAttractors() was not able to identify any attractors! Please check the supplied parameters and restart!")
    
    if (length(network$genes) %% 32 == 0)
      numElementsPerEntry <- as.integer(length(network$genes) / 32)
    else
      numElementsPerEntry <- as.integer(length(network$genes) / 32  + 1)
    
    if (!is.null(result$stateInfo))
    {
      result$stateInfo$table <- matrix(result$stateInfo$table,nrow=numElementsPerEntry)
    
      if (!is.null(result$stateInfo$initialStates))
        result$stateInfo$initialStates <- matrix(result$stateInfo$initialStates,nrow=numElementsPerEntry)
    }
    
    for (i in seq_len(length(result$attractors)))
    {
      result$attractors[[i]]$involvedStates <- matrix(result$attractors[[i]]$involvedStates,nrow=numElementsPerEntry)
      if (canonical)
      # reorder states
        result$attractors[[i]]$involvedStates <- canonicalStateOrder(result$attractors[[i]]$involvedStates)
        
      if (!is.null(result$attractors[[i]]$initialStates))
        result$attractors[[i]]$initialStates <- matrix(result$attractors[[i]]$initialStates,nrow=numElementsPerEntry)
      if (!is.null(result$attractors[[i]]$nextStates))
        result$attractors[[i]]$nextStates <- matrix(result$attractors[[i]]$nextStates,nrow=numElementsPerEntry)
        
      if (result$attractors[[i]]$basinSize == 0)
        result$attractors[[i]]$basinSize <- NA
    }
    
    # order attractors according to their lengths
    attractorLengths <- sapply(result$attractors,function(attractor)ncol(attractor$involvedStates))  
    reordering <- order(attractorLengths)
    result$attractors <- result$attractors[reordering]
    
    if (!is.null(result$stateInfo))
    {
      inverseOrder <- sapply(seq_along(reordering),function(x)which(reordering == x))    
      result$stateInfo$attractorAssignment <- inverseOrder[result$stateInfo$attractorAssignment]
    }
    
    # extend the resulting structure by additional information, and assign a class      
    result$stateInfo$genes <- network$genes
    result$stateInfo$fixedGenes <- network$fixed

    if (!is.null(result$stateInfo$table))
      class(result$stateInfo) <- "BooleanStateInfo"
    class(result) <- "AttractorInfo"
    return(result)
  }
}
