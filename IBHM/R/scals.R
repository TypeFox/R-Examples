# Scalarization functions ##############################

ScalFunctions <- function(filter = NULL){
  # Dot product scalarization
  DotPr <- function(x,d){  
    return(x%*%d[-1]+d[[1]])
  }
  expr(DotPr) <- 'dot.pr'



  # Radial scalarization
  Radial <- function(x,d){
    n <- dim(x)[[2]]
    m <- dim(x)[[1]]
    
    d_m <- matrix(rep(d[-1],m),nrow=m,ncol=n)
    
    rowSums(d[[1]]*(x-d_m)^2)
  }
  expr(Radial) <- 'radial'

  # Radial scalarization - euclidean distance
  RootRadial <- function(x,d){
    n <- dim(x)[[2]]
    m <- dim(x)[[1]]
    
    d_m <- matrix(rep(d[-1],m),nrow=m,ncol=n)
    
    rowSums(d[[1]]*sqrt((x-d_m)^2))
  }
  expr(RootRadial) <- 'root.radial'
  
  res <- list(DotPr, Radial, RootRadial)  
  if(!is.null(filter)){
    res <- Filter(function(v) any(expr(v)==filter), res)
  }
  if(length(res)==0){
    stop('Invalid scal list given : ', paste(filter,sep=', '))
  }
  
  res
}

ScalParamDim <- function(scal, x){
  switch(expr(scal), 
         dot.pr = ncol(x)+1,
         radial = ncol(x)+1,
         root.radial = ncol(x)+1,
         stop('Unknown scal: ',expr(scal)))
}