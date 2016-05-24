
# Convert a SymbolicBooleanNetwork object 
# <network> to a BooleanNetwork object
symbolicToTruthTable <- function(network)
{
  stopifnot(inherits(network, "SymbolicBooleanNetwork"))
  res <- list(genes = network$genes,
              fixed = network$fixed,
              interactions = lapply(network$interactions, function(int)
              {
                newInt <- .Call("getTruthTable_R", int, length(network$genes))
                names(newInt) <- c("input","func")
                newInt$expression <- stringFromParseTree(int)
                return(newInt)
              }))
  class(res) <- "BooleanNetwork"
  return(res)
}
