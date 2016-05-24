create.parameter.df <-
function(alpha,x,y,e)
  {
   par_output <- c(alpha, x, y, e)
   names(par_output) <- c("alpha", "x", "y", "e")
   return(as.data.frame(par_output))
  }
