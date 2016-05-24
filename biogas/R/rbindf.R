# Modified: 22 JUNE 2015
# From monitoR package

rbindf <- function(...) {
  
  l <- list(...)
  if(length(l)==1) l <- l[[1]]
  nn <- length(l)

  x <- l[[1]]
  if(length(l)>1){
      for(i in 2:nn) {
        y <- l[[i]]
        if(!all(yinx <- names(y) %in% names(x))) {
          x[, names(y)[!yinx]] <- NA
        } 
        if(!all(xiny <- names(x) %in% names(y))) {
           y[, names(x)[!xiny]] <- NA
        } 
        x <- rbind(x, y)
      }
  }
  return(x)
}


