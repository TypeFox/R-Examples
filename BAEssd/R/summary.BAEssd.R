summary.BAEssd <-
function(object,...){
  ### Obtain the relevant information
  out <- as.numeric(object$history[object$history$n==object$n,])
  
  names(out) <- colnames(object$history)
  attributes(out) <- c(attributes(out),attributes(object$n))
  class(out) <- "summary.BAEssd"

  return(out)
}
