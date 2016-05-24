`setplotmat` <-
function(x,y)
  {
    ###  given an x and y vectors of lat-lon (e.g.)
    ###  return a matrix of points
    EX =  matrix( rep(x, length(y)), byrow=TRUE, ncol=length(x))
    WHY =  matrix( rev(rep(y, length(x))), byrow=FALSE, ncol=length(x))
    
    return(list(x=EX, y=WHY))
  }

