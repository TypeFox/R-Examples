"pp" <-
function(x, A=NULL, B=NULL, a=0, sort=TRUE) {

   if(! is.null(a)) {
      if(a < 0 | a > 0.50) {
         warning("Plotting position parameter a is invalid, not in [0,0.5]")
         return()
      }
      A <- -a
      B <- 1 - 2*a
   }
   if(is.null(A)) {
      warnings("Plotting position parameter A is NULL")
      return(NULL)
   }
   if(is.null(B)) {
      warnings("Plotting position parameter B is NULL")
      return(NULL)
   }
   if(A <= -1 | A >= B) {
      warnings("Plotting position parameters A or B are invalid")
      return(NULL)
   }
    
   denom <- length(x) + B
   ranks <- rank(x, ties.method="first")
    
   if(sort) {
      return( (sort(ranks) + A) / denom)
   } else {
      return(      (ranks  + A) / denom)
   }
}

