
dscore.default <- function(x, ...)
 {
   if(any(x>1 | x<0))
    stop("argument 'x' should be a vector containing probabilities")
   ans <- .Call("scorefreqs", x, PACKAGE="SNPassoc")
   names(ans) <- paste(0:(2*length(x)), "alleles")
   ans 
 }


