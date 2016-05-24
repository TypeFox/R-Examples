`ojaMedianEx` <- function(X, control=ojaMedianControl(...), ...){

    action <- 1
    param2 <- param3 <- param4 <- debug <- 0
    rows <- dim(X)[1]
    cols <- dim(X)[2]
    outvec <- c(1:cols)
    x <- y <- 1
   storage.mode(rows) <- "integer"
   storage.mode(cols) <- "integer"
   storage.mode(X) <- "double"
   storage.mode(outvec) <- "double"
   
    res<-.C("r_oja", rows, cols, X, vec = outvec, y, as.integer(action), as.double(control$maxlines), as.double(param2), as.integer(param3), as.integer(param4), as.integer(debug),1)
    RES <- res$vec
  names(RES)<-colnames(X)
  return(RES)
  }
