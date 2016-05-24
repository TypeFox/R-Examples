# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
summary.cont.cont <- function(object, ...){

  cat("\nThe forecasts are continuous, the observations are continous.\n")
if(object$baseline.tf){cat("Baseline data provided. \n")} else{
  cat("Sample baseline calcluated from observations.\n")}

  cat(paste("MAE               = ", formatC(object$MAE, digits = 4), "\n"))
  cat(paste("ME                = ", formatC(object$ME, digits = 4), "\n"))
  cat(paste("MSE               = ", formatC(object$MSE, digits = 4), "\n"))
  cat(paste("MSE - baseline    = ", formatC(object$MSE.baseline, digits = 4), "\n"))
  cat(paste("MSE - persistence = ", formatC(object$MSE.pers, digits = 4), "\n"))
  cat(paste("SS  - baseline     = ", formatC(object$SS.baseline, digits = 4), "\n"))

}


