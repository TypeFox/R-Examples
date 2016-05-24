dectobin <-
function(x,forward=TRUE,nl=NULL){
  if(forward){
    xi <- as.integer(x)
    N <- length(x)
    xMax <- max(x)	
    ndigits <- (floor(logb(xMax, base=2))+1)
    Base.b <- array(NA, dim=c(N, ndigits))
    for(i in 1:ndigits){#i <- 1
      Base.b[, ndigits-i+1] <- (x %% 2)
      x <- (x %/% 2)
    }
    tmp <- if(N ==1) rev(Base.b[1, ]) else rev(Base.b)
    if(length(tmp)<nl) tmp <- c(tmp,rep(0,nl-length(tmp)))
    return(tmp)
  }else{
    sum(2^(which((x))-1))
  }
}

