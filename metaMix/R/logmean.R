logmean <- function(logL) {


  my.max <- max(logL)
  my.mod.logL <- logL - my.max
  return( log(mean(exp(my.mod.logL))) + my.max)
  
  
}
