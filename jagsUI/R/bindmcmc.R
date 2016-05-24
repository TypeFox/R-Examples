
bind.mcmc <- function(mcmc.list1,mcmc.list2,start,n.new.iter){
  
  nchains <- length(mcmc.list1)
  
  samples <- list()
  
  for (i in 1:nchains){
    
    d <- rbind(mcmc.list1[[i]],mcmc.list2[[i]])
    
    samples[[i]] <- mcmc(data=d,start=start,end=(end(mcmc.list1[[i]])+n.new.iter),thin=thin(mcmc.list1[i]))
    
  }
  
  return(as.mcmc.list(samples))
  
  
}