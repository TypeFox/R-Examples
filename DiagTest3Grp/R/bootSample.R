########################################################################
#This function draws a randomly sampled vector based on inputted vector x
#########################################################################
bootSample <-
function(x,seed0)  
  { 
    n <- length(x)
    if(!missing(seed0) &!is.null(seed0)) set.seed(seed0)
    idx <- sample(1:n,n,replace=T)
    new.x <- x[idx]
    return(new.x)
  }

