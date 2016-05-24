makeCor <-
function(x) {
 if (is.matrix(x)) stop("x is not a vector.")
 upper       <- matrix(x,ncol=10, byrow=FALSE)
 diag(upper) <- 0
 lower <- matrix(x,ncol=10, byrow=TRUE)
 res <- lower + upper
 return(res)
 }