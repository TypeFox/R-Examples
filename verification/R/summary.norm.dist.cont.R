# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
summary.norm.dist.cont <- function(object, ...){

  cat("\nThe forecasts are a normal probability distribution. \n")
  cat("The observations are continuous.\n\n")
  cat(paste("Average crps score      = ", formatC(object$CRPS, digits = 4), "\n"))
  cat(paste("Average ignorance score = ", formatC(object$IGN, digits = 4), "\n"))
}
