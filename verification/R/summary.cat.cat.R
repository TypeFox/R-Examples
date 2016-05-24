# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
summary.cat.cat <- function(object, ...){

  cat("\nThe forecasts are categorical, the observations are categorical.\n")

  cat(paste("Percent Correct      = ", formatC(object$pc, digits = 2), "\n"))
  cat(paste("Heidke Skill Score   = ", formatC(object$hss, digits = 3), "\n"))
  cat(paste("Pierce Skill Score   = ", formatC(object$pss, digits = 3), "\n"))
  cat(paste("Gerrity Score        = ", formatC(object$gs, digits = 3), "\n"))
  cat("\n Statistics considering each category in turn. \n \n" )
  cat(c("Threat Score             ", formatC(object$ts, digits = 3), "\n"))
  cat(c("Bias by cat.             ",   formatC(object$bias2, digits = 3), "\n"))
  cat(c("Percent correct by cat.  ",   formatC(object$pc2, digits = 3), "\n") )
  cat(c("Hit Rate (POD) by cat.   ",   formatC(object$h, digits = 3), "\n") )
  cat(c("False Alarm Rate by cat. ",   formatC(object$f, digits = 3), "\n") )
  cat(c("False Alarm Ratio by cat.",   formatC(object$false.alarm.ratio, digits = 3), "\n") )
#  
}
