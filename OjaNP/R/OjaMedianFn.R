`ojaMedianFn` <-
function(X,x)
    { y <- 1
      param1 <- param2 <- param3 <- param4 <- debug <- 0
      rows <- dim(X)[1]
      cols <- dim(X)[2]
      outvec <- c(1:cols)
      storage.mode(rows) <- "integer"
      storage.mode(cols) <- "integer"
      storage.mode(X) <- "double"
      storage.mode(outvec) <- "double"
      action <- 3
      storage.mode(x) <- "double"
      res<-.C("r_oja",rows,cols,X,vec=outvec,y,as.integer(action),as.double(x),as.double(param2),as.integer(param3),as.integer(param4),as.integer(debug))
      (res$vec)[1]
    }
