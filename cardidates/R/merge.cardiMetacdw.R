merge.cardiMetacdw <-
function(x, y, ...) {
  if (!is(y, "cardiMetacdw")) 
    stop("incompatible argument. 'y' must be of class 'cardiMetacdw'")
  if (nrow(y$metares) != 1)
    stop("sorry, second argument must contain fit with only one sample")
  
  ndx <- which(as.character(x$metares$sample) == 
               as.character(y$metares$sample))
  if (length(ndx) == 0) {
    ## add new fit to the objec
    x$metares <- rbind(x$metares, y$metares)
    j <- length(x$weibullfits)
    x$weibullfits[j+1] <- y$weibullfits
  } else {
    ## replace existing fit
    x$metares[ndx,]    <- y$metares
    x$weibullfits[ndx] <- y$weibullfits
  }
  #class(x) <- c("list", "cardiMetacdw")
  return(x)
}

