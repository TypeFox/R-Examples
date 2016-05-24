rescale <- 
  function(x, ...){
    
    if(is.null(attr(x, "class"))){
      return(x)
    } 
    else  UseMethod("rescale", x)
    
  }