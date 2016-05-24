`ojaMedianGrid` <-
function(X, control=ojaMedianControl(...), ...){
    
    action <- 2
    param4 <- debug <- 0

   rows <- dim(X)[1]
   cols <- dim(X)[2]
   outvec <- c(1:cols)
   SEED <- sample(1:5000)
   x <- y <- 1
   storage.mode(rows) <- "integer"
   storage.mode(cols) <- "integer"
   storage.mode(X) <- "double"
   storage.mode(outvec) <- "double"
   
    solution<-.C("r_oja", rows, cols, X, vec = outvec, y, as.integer(action), as.double(control$eps), as.double(control$chi2), as.integer(control$samples), as.integer(param4), as.integer(debug))
    #output <- solution$vec
    output <- solution$vec
    return(output)
  }
