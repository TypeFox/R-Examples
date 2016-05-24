stateTransition <- function(network,state,type=c("synchronous","asynchronous","probabilistic"),
                            geneProbabilities, chosenGene, chosenFunctions, timeStep=0)
{
  stopifnot(inherits(network,"BooleanNetwork") | inherits(network,"SymbolicBooleanNetwork") 
            | inherits(network,"ProbabilisticBooleanNetwork"))

  type <- match.arg(type, c("synchronous","asynchronous","probabilistic"))
  if (inherits(network,"SymbolicBooleanNetwork"))
  {
    if (type != "synchronous")
      stop("Only synchronous updates are allowed for symbolic networks!")
      
    if (!is.null(dim(state)))
    {
      if (ncol(state) != length(network$genes))
        stop(paste("\"state\" must be either a vector with one value for each gene,",
                   "or a matrix with the genes in the columns and multiple predecessor states in the rows!"))
      state <- as.integer(as.matrix(state[nrow(state):1,]))
    }
    else
    {
      if (!is.numeric(state) || length(state) != length(network$genes))
        stop(paste("\"state\" must be either a vector with one value for each gene,",
                   "or a matrix with the genes in the columns and multiple predecessor states in the rows!"))
      state <- as.integer(state)
    }

    if (is.null(network$internalStructs) || checkNullPointer(network$internalStructs))
    # refresh internal tree representations if necessary
      network$internalStructs = .Call("constructNetworkTrees_R",network);
      
    on.exit(.C("freeAllMemory", PACKAGE = "BoolNet"))  
    res <- .Call("symbolicStateTransition_R",
                 network$internalStructs, 
                 state,
                 as.integer(timeStep))
  }
  else
  {
    if (length(state) != length(network$genes))
        stop("The state must consist of exactly one value for each gene!")
            
    nonFixedIndices = (network$fixed == -1)

    res <- state
    
    if (inherits(network,"BooleanNetwork") & type ==  "probabilistic")
      type <- "synchronous"
      
    if (!inherits(network,"BooleanNetwork") & length(type) == 3)
      type <- "probabilistic"
    
    if (type ==  "probabilistic")
    {
      if (missing(chosenFunctions) || is.null(chosenFunctions))
      {
        chosenFunctions <- sapply(network$interactions,function(gene)
                                  {
                                    distr <- c(0,cumsum(sapply(gene,function(func)func$probability)))
                                    r <- runif(n=1)
                                    idx <- 0
                                    for (i in seq_along(gene))
                                    {
                                      if (r > distr[i] & r <= distr[i+1])
                                      {
                                        idx <- i
                                        break
                                      }
                                    }
                                    idx
                                  })
      }
      else
        if (length(chosenFunctions) != length(network$genes))
          stop("Please provide a function index for each gene!")    

      res[nonFixedIndices] <- mapply(function(i,f)
                { 
                  if(network$interactions[[i]][[f]]$input[1] == 0)
                  # this is a constant gene with no transition function
                    return(network$interactions[[i]][[f]]$func[1])
                  input = state[network$interactions[[i]][[f]]$input]
                  return(network$interactions[[i]][[f]]$func[bin2dec(rev(input),length(input)) + 1])
                }, which(nonFixedIndices), chosenFunctions[nonFixedIndices])                                        
    }
    else
    {
      if (!inherits(network,"BooleanNetwork"))
        stop("Please choose type=\"probabilistic\" for a probabilistic Boolean network!")
      changeIndices <- switch(type,
        synchronous = which(nonFixedIndices),
        asynchronous =
        {
          if (missing(chosenGene) || is.null(chosenGene))
          {
            if (missing(geneProbabilities) || is.null(geneProbabilities))
              sample(which(nonFixedIndices),1)
            else
            {
              if (length(geneProbabilities) != length(network$genes))
                stop("Please supply exactly one probability for each gene!")
                
              if (abs(1.0 - sum(geneProbabilities)) > 0.0001)
                stop("The supplied gene probabilities do not sum up to 1!")
                
              
              if (sum(geneProbabilities[nonFixedIndices]) < 0.0001)
              	stop("There is no non-fixed gene with a probability greater than 0!")  
                
              geneProbabilities[nonFixedIndices] <- geneProbabilities[nonFixedIndices]/sum(geneProbabilities[nonFixedIndices])
              
              if (sum(nonFixedIndices) != length(network$genes))
                geneProbabilities[!nonFixedIndices] <- 0 
              distr <- c(0,cumsum(geneProbabilities))
              r <- runif(n=1)
              idx <- 0
              for (i in seq_along(network$genes))
              {
                if (r > distr[i] & r <= distr[i+1])
                {
                  idx <- i
                  break
                }
              }
              idx
            }
          }
          else
          {
            if (is.character(chosenGene))
            {
              chosenGeneIdx <- which(network$genes == chosenGene)
              if (length(chosenGeneIdx) == 0)
                stop(paste("Gene \"",chosenGene,"\" does not exist in the network!",sep=""))
              chosenGeneIdx
            }
            else
              chosenGene
          }
      })

      res[changeIndices] <- sapply(changeIndices,function(i)
      {    
        if(network$interactions[[i]]$input[1] == 0)
        # this is a constant gene with no transition function
          return(network$interactions[[i]]$func[1])
        input = state[network$interactions[[i]]$input]
        return(network$interactions[[i]]$func[bin2dec(rev(input),length(input)) + 1])
      })
    }

    res[!nonFixedIndices] = network$fixed[!nonFixedIndices]
  }
  
  names(res) <- network$genes
  return(res)
}
