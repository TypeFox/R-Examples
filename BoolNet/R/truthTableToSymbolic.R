
# Convert a BooleanNetwork object 
# <network> to a SymbolicBooleanNetwork object
truthTableToSymbolic <- function(network, generateDNFs=FALSE)
{
  stopifnot(inherits(network, "BooleanNetwork"))
  
  res <- list(genes = network$genes,
              fixed = as.integer(network$fixed),
              interactions = lapply(network$interactions, function(int)
              {
                func <- tryCatch(
                {
                  if (generateDNFs != FALSE)
                    parseBooleanFunction(getDNF(int$func, network$genes[int$input], mode=generateDNFs),
                                              network$genes)
                  else                                              
                    parseBooleanFunction(int$expression,
                                              network$genes)
                },
                error =  function(e)
                {
                  parseBooleanFunction(getDNF(int$func, network$genes[int$input], mode=generateDNFs),
                                            network$genes)
                })
                return(func)
              }))

  names(res$fixed) <- res$genes
  res$internalStructs <- .Call("constructNetworkTrees_R",res)
  res$timeDelays <- apply(sapply(res$interactions,maxTimeDelay,genes=res$genes),1,max)
              
  class(res) <- "SymbolicBooleanNetwork"
  return(res)
}
