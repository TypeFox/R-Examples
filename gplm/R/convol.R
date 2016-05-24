"convol" <- function(x,h=1,grid=NULL,y=1,w=1,p=2,q=2,product=TRUE,sort=TRUE){

  x  <- as.matrix(x)        ## nxd
  h  <- as.matrix(h)        ## scalar / 1xd ( future: nx1 / nxd )
  w  <- as.matrix(w)        ## scalar / nx1 / 1xm / nxm
  y  <- as.matrix(y)        ## nxc
  
  n <- nrow(x)
  d <- ncol(x)
  c <- ncol(y)
  
  if(!(nrow(y) %in% c(1,n))) stop("convol: y must have 1 or nrow(x) rows")
  if (nrow(y)==1){ y <- t(t(matrix(1, n, c))*as.vector(y)) }
  
  if (sort){
    or <- order(x[,1])
    ro <- order((1:n)[or])
    x  <- x[or,,drop=FALSE]
    y  <- y[or,,drop=FALSE]
    
    if (nrow(h)==nrow(x)){ h <- h[or,,drop=FALSE] }
    if (nrow(w)==nrow(x)){ w <- w[or,,drop=FALSE] }
  }
  
  if (is.null(grid)){
    grid <- x
    havegrid <- FALSE
  }else{
    grid <- as.matrix(grid)     ## mxd
    havegrid <- TRUE
  }
  
  m <- nrow(grid)
  r <- matrix(NA,m,c)
  
  if (havegrid && sort){
    or.grid <- order(grid[,1])
    ro.grid <- order((1:m)[or.grid])
    grid  <- grid[or.grid,,drop=FALSE]
    
    if (ncol(w)==nrow(grid)){ w <- w[,or.grid,drop=FALSE] }
  }
  
  ## h: scalar or 1xd
  if (nrow(h)==1){
    if (ncol(h)!=d){ h <- matrix(h,1,d) }
  }

  ##print("h: scalar or 1xd")
  x  <- t( t(x) /as.vector(h) ) 
  grid <- t( t(grid)/as.vector(h) )
  
  if (ncol(w)==1){              ## w: scalar or nx1
    ##print("w: scalar or nx1")
    
    dim <- as.double(c(n,m,d,c,p,q,product+0))
    r  <- .C("convol",dim,as.double(x),as.double(y*as.vector(w)),as.double(grid),as.double(rep(0,m*c)),PACKAGE="gplm")
    r  <- matrix(r[[5]],nrow=m)/prod(h)
  }else{                        ## w: 1xm or nxm
    ##print("w: 1xm or nxm")
    if (nrow(w)==1){
      dim <- as.double(c(n,m,d,c,p,q,product+0))
      r  <- .C("convol",dim,as.double(x),as.double(y),as.double(grid),as.double(rep(0,m*c)),PACKAGE="gplm")
      r  <- matrix(r[[5]],nrow=m)/prod(h)
      r  <- t( t(r)*as.vector(w) )
    }else{
      for (j in 1:m){
        dim <- as.double(c(n,1,d,c,p,q,product+0))
        rj <- .C("convol",dim,as.double(x),as.double(w[,j]*y),as.double(grid[j,]),as.double(rep(0,c)),PACKAGE="gplm")
        r[j,] <- rj[[5]]/prod(h)
      }
    }
    
  }

  return(r)
}
  
  
