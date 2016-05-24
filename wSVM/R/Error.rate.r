Error.rate = function(Fit, Y){
   if(is.vector(Y)) N <- length(Y)
   else N <- nrow(Y)
   Error = sum(Fit != Y) / N
   return(Error)
}
