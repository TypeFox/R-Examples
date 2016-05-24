# Creates a deterministic Boolean network from a Probabilistic Boolean network
# by choosing the functions specified in <functionIndices> from the
# network <probabilisticNetwork>.
chooseNetwork <- function(probabilisticNetwork,functionIndices,dontCareValues=NULL,readableFunctions=FALSE)
{
  stopifnot(inherits(probabilisticNetwork,"ProbabilisticBooleanNetwork") 
        | inherits(probabilisticNetwork,"BooleanNetworkCollection"))

  if (length(functionIndices) != length(probabilisticNetwork$genes))
    stop("Please provide a vector of function indices for each gene!")

  if (inherits(probabilisticNetwork,"ProbabilisticBooleanNetwork"))
  {
    interactions <- mapply(function(interaction,index)
          {
            list(input=interaction[[index]]$input,
              func=interaction[[index]]$func,
              expression=interaction[[index]]$expression)
          },
          probabilisticNetwork$interactions,functionIndices,SIMPLIFY=FALSE)
  }
  else
  {
    if (!is.null(dontCareValues) && 
        length(dontCareValues) == length(probabilisticNetwork$genes) &&
        is.null(names(dontCareValues)))
          names(dontCareValues) <- probabilisticNetwork$genes
    interactions <- mapply(function(interaction,index,gene)
          {
            func <- interaction[[index]]$func
            dcPos <- which(func == -1)
            
            if (length(dcPos) > 0)
            {
              if (is.null(dontCareValues) || is.null(dontCareValues[[gene]]))
                stop(paste("No values for the \"don't care\" ",
                           "entries were specified for gene \"",
                           gene,"\"!",sep=""))

              if (!all(dontCareValues[[gene]] %in% c(0,1)))
                stop(paste("Invalid values for \"don't care\" entries specified for gene \"",
                           gene,"\"!",sep=""))
                           
              if (length(dontCareValues[[gene]]) != length(dcPos))
                stop(paste("There must be exactly ",length(dcPos),
                           " value(s) for \"don't care\" entries in the function for gene \"",
                           gene,"\"!",sep=""))

              func[[dcPos]] <- dontCareValues[[gene]]
              expression <- getInteractionString(readableFunctions,
                                                 func,
                                                 probabilisticNetwork$genes[interaction[[index]]$input])
            }
            else
              expression <- interaction[[index]]$expression

            list(input=interaction[[index]]$input,
              func=func,
              expression=expression)
          },
          probabilisticNetwork$interactions,
          functionIndices,
          probabilisticNetwork$genes,SIMPLIFY=FALSE)

  }
  
  res <- list(genes=probabilisticNetwork$genes,
      interactions=interactions,
      fixed=probabilisticNetwork$fixed)
  class(res) <- "BooleanNetwork"
  return(res)
}
