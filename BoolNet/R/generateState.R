generateState <- function(network, specs, default=0)
{
  stopifnot(inherits(network,"ProbabilisticBooleanNetwork")
            || inherits(network,"BooleanNetwork")
            || inherits(network,"SymbolicBooleanNetwork"))
  
  if (!all(names(specs) %in% network$genes))
    stop(paste("Undefined gene names:",
               paste(setdiff(names(specs), network$genes), collapse=", ")))
  
  if (!all(unlist(specs) %in% c(0,1)))
    stop("Please provide only Boolean values!")
 
 lengths <- unique(sapply(specs, length))
 if (length(lengths) > 1) 
    stop("The number of specifications for each gene must be the same!")
 
 if (lengths == 1)
 {   
   state <- rep(default, length(network$genes))
   names(state) <- network$genes
   state[names(specs)] <- specs
 }
 else
 {
   state <- matrix(rep(default, length(network$genes) * lengths), nrow=lengths)
   colnames(state) <- network$genes
   
   for (i in seq_along(specs))
    state[,names(specs)[i]] <- specs[[i]]
 }
 return(state)                
}
