## Convert a rbugs object to a coda object
rbugs2coda <- function(model, burnin = NULL, thin = NULL)
{
    
    if( class(model) != "rbugs") stop("\nThe model is not a rbugs obejct.\n")
    chain.name <- c()
    for (i in 1:model$n.chain)
      chain.name <- c(chain.name, paste("chain",i,sep=""))

    obj <- list()
    for (i in 1:model$n.chain)
      obj <- c(obj,list(as.mcmc(model[[chain.name[i]]])))
    
    if(length(burnin) != 0){
      for (i in 1:model$n.chain)
       obj[[i]] <- mcmc(obj[[i]], start = burnin, end = nrow(obj[[i]])) 
    } 
    
    if(length(thin) != 0){
      for (i in 1:model$n.chain)
       obj[[i]] <- mcmc(obj[[i]], start = 1, end = nrow(obj[[i]]), thin = thin) 
    }
    
    obj <- as.mcmc.list(obj)
    
    obj
}
