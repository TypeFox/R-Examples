# Map (x,y) to (newx,newy)
# x    : a vector of 1 * n
# y    : a vector of 1 * n or a n * p matrix
# newx : vector of 1 * m
# newy : vector of 1 * m or a matrix of m * optns$
  
  # if( is.vector(y) ){
    # return(y[is.element(x,newx)])
  # }else if(is.matrix(y)){
    # return(y[is.element(x,newx),])
  # }else{
    # warning('y cannot be empty!\n')
    # return(NaN)
  # }   
# }

MapX1D <- function(x, y, newx) {
    # if (!all(newx %in% x)) 
        # warning('Interpolation occured: you might want to increase the obsGrid coverage')
        
    # if (min(newx) + 100 * .Machine$double.eps < min(x) || max(newx) > max(x) + 100 * .Machine$double.eps)
       # warning('Extrapolation occured')
    if (is.vector(y)){
        # newy <- approxExtrap(x, y, newx, method='linear')$y
      newy <- approx(x, y, newx, method='linear')$y
    } else {
        # newy <- apply(y, 2, function(yy) approxExtrap(x, yy, newx, method='linear')$y)
      newy <- apply(y, 2, function(yy) approx(x, yy, newx, method='linear')$y) 
    }
    if (any(is.nan(newy))){
      stop('NA \'s during the mapping from(x,y) to (newx,newy)')
    }
    
    return(newy)
}
