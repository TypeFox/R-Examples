
### The function tests whether a given gradient is given
### observation-wise.  It tests essentially the # of rows
### in the gradient
observationGradient <- function(g, nParam) {
   if(is.null(dim(g))) {
      if(nParam == 1 & length(g) > 1)
          return(TRUE)
      return(FALSE)
   }
   if(nrow(g) == 1)
       return(FALSE)
   return(TRUE)
}
