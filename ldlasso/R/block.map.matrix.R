block.map.matrix <-
function( block.cood ){
  X <- matrix(FALSE, nrow = length(block.cood)-1, ncol = length(block.cood)-1 )
  i <- 1
  while( i < length(block.cood) ){  
    j <- i+1
    if( block.cood[j]==1 ){
      # do nothing
    }else{
      while( block.cood[j] == 0 ){
        X[i,j] <- TRUE
        j <- j+1
      }
    }
    i <- i+1
  }
  return(t(X))
}

