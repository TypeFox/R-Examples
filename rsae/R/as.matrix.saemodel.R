as.matrix.saemodel <-
function(x, ...){
   areaID <- x$areaID
   X <- x$X
   y <- x$y
   res <- cbind(y, X, areaID)
   res <- as.matrix(res)
   return(res)
}

