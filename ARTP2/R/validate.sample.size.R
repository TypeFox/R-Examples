
validate.sample.size <- function(family, lambda, ncases, ncontrols, nsamples){
  
  if(!is.list(nsamples)){
    msg <- 'nsamples should be a list'
    stop(msg)
  }
  
  if(!is.list(ncases)){
    msg <- 'ncases should be a list'
    stop(msg)
  }
  
  if(!is.list(ncontrols)){
    msg <- 'ncontrols should be a list'
    stop(msg)
  }
  
  if(family == 'gaussian'){
    
    if(length(nsamples) == 0){
      msg <- 'nsamples cannot be an empty list when when family is \'gaussian\''
      stop(msg)
    }
    
    if(length(ncases) > 0){
      msg <- 'Users should not specify ncases when family is \'gaussian\''
      stop(msg)
    }
    
    if(length(ncontrols) > 0){
      msg <- 'Users should not specify ncontrols when family is \'gaussian\''
      stop(msg)
    }
    
    if(length(nsamples) != length(lambda)){
      msg <- 'Length of nsamples and lambda should be equal'
      stop(msg)
    }
    
  }else{ # family == 'binomial
    
    if(length(nsamples) > 0){
      msg <- 'Users should not specify nsamples when family is \'binomial\''
      stop(msg)
    }
    
    if(length(ncases) == 0){
      msg <- 'ncases cannot be an empty list when when family is \'binomial\''
      stop(msg)
    }
    
    if(length(ncontrols) == 0){
      msg <- 'ncontrols cannot be an empty list when when family is \'binomial\''
      stop(msg)
    }
    
    if(length(ncases) != length(lambda)){
      msg <- 'Length of ncases and lambda should be equal'
      stop(msg)
    }
    
    if(length(ncontrols) != length(lambda)){
      msg <- 'Length of ncontrols and lambda should be equal'
      stop(msg)
    }
    
    for(i in 1:length(ncases)){
      if(length(ncases[[i]]) != length(ncontrols[[i]])){
        msg <- 'Length of each element of ncases and ncontrols should be equal'
        stop(msg)
      }
    }
  }
  
}
