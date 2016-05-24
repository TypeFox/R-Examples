"normalize.vector" <-
  function(x, mass=NULL) {
    x <- as.matrix(x);
    dx <- dim(x); 
    
    if(dx[2]==1)
      x <- as.numeric(x)
    
    if(!is.null(mass)) {
      if (dx[1] != (length(mass)*3))
        stop("normalize.vector: incorrect length of mass")
    }

    if(is.matrix(x))
      return(t( t(x) / sqrt(inner.prod(x,x,mass))))
    else
      return(x / sqrt(inner.prod(x,x,mass))) 
  }
