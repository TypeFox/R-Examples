
gen.inits <- function(inits,n.chains){
  
  if(is.null(inits)){
    init.values <- vector("list",length=n.chains)
    for(i in 1:n.chains){     
      init.values[[i]]$.RNG.name="base::Mersenne-Twister"
      init.values[[i]]$.RNG.seed=abs(.Random.seed[i+1])
    }
  }else if(is.list(inits)){
    if(length(inits)!=n.chains){stop('Length of initial values list != number of chains')}
    init.values <- inits
    for(i in 1:n.chains){
      init.values[[i]]$.RNG.name="base::Mersenne-Twister"
      init.values[[i]]$.RNG.seed=abs(.Random.seed[i+1])
    }
    
  } else if (is.function(inits)){
    init.values <- list()
    for (i in 1:n.chains){
      init.values[[i]] <- inits()
      init.values[[i]]$.RNG.name="base::Mersenne-Twister"
      init.values[[i]]$.RNG.seed=abs(.Random.seed[i+1])
    }
    
  } else {stop('Invalid initial values. Must be a function or a list with length=n.chains')}
    
  return(init.values) 
  }
