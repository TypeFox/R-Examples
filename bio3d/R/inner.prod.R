"inner.prod" <-
  function(x,y,mass=NULL) {
    x <- as.matrix(x); y <- as.matrix(y);
    dx <- dim(x); dy <- dim(y);
    
    if(dx[1]!=dy[1])
      stop("inner.prod: unequal vector lengths")
    
    if(dx[2]>1 && dy[2]>1) {
      if(dx[2]!=dy[2])
        stop("inner.prod: unequal vector lengths")
    }
    
    if(dx[2]==1)
      x <- as.numeric(x)
    if(dy[2]==1)
      y <- as.numeric(y)
    
    if(!is.null(mass)) {
      if (dx[1] != (length(mass)*3)) 
        stop("inner.prod: incorrect length of mass")
    }
    
    if(is.null(mass))
      mass <- 1
    else
      mass <- rep(mass,each=3)
    
    if(is.matrix(x) || is.matrix(y))
      return(colSums((x*y)*mass^2))
    else
      return(sum(x*y*mass^2))
  }
