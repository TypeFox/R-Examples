print.incanc <- function(x, ...){
# x 
  cond1 <- !is.null(x$noise[[2]])
  if (cond1){
    cond2 <- length(x$noise[[2]])==0
    if (cond2){
       cat("-- No noise units considered -- \n")
    }else{
       cat("Noise units: \n", x$noise[[2]], "\n")
    }
  }
 
 
   cat ("\n---INCA index to estimate the number of clusters---\n")
   cat ("Clustering method: ", x$method, " \n")                               
   
 
   x <- x$INCAindex                                 
   k <- length(x)
   for (i in 2:k){
      if(!is.na(x[i])){ 
         cat(" k= ", i,"      ",  format(x[i], digits=2), " \n", sep="")
      }
   }
   
   
}
