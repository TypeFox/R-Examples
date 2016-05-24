# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
summary.prob.bin <- function(object, ...){

  cat("\nThe forecasts are probabilistic, the observations are binary.\n")
if(object$baseline.tf){cat("Baseline data provided. \n")} else{
  cat("Sample baseline calculated from observations.\n")}
  cat(paste("Brier Score (BS)           = ", formatC(object$bs, digits = 4), "\n"))
  cat(paste("Brier Score - Baseline     = ", formatC(object$bs.baseline, digits = 4), "\n"))
  cat(paste("Skill Score                = ", formatC(object$ss, digits = 4), "\n"))
  cat(paste("Reliability                = ", formatC(object$bs.reliability, digits = 4), "\n"))
  cat(paste("Resolution                 = ", formatC(object$bs.resol, digits = 4), "\n"))
  cat(paste("Uncertainty              = ", formatC(object$bs.uncert, digits = 4), "\n")) 
}
