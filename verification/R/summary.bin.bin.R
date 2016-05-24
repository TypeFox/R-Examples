# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
summary.bin.bin <- function(object, ...){
  ## print function for binary forecast, binary outcome
  cat("\nThe forecasts are binary, the observations are binary.\n")
  cat("The contingency table for the forecast \n")
  print(object$tab)
  cat("\n")
 
  cat(paste("PODy = ", formatC(object$POD, digits = 4), "\n"))
  cat(paste("Std. Err. for POD = ", formatC(object$POD.se), "\n"))
  cat(paste("TS   = ", formatC(object$TS, digits = 4), "\n"))
  cat(paste("Std. Err. for TS = ", formatC(object$TS.se), "\n"))
  cat(paste("ETS  = ", formatC(object$ETS, digits = 4), "\n"))
  cat(paste("Std. Err. for ETS = ", formatC(object$ETS.se), "\n"))
  cat(paste("FAR  = ", formatC(object$FAR, digits = 4), "\n"))
  cat(paste("Std. Err. for FAR = ", formatC(object$FAR.se), "\n"))
  cat(paste("HSS  = ", formatC(object$HSS, digits = 4), "\n"))
  cat(paste("Std. Err. for HSS = ", formatC(object$HSS.se), "\n"))
  cat(paste("PC   = ", formatC(object$PC, digits = 4), "\n"))
  cat(paste("Std. Err. for PC = ", formatC(object$PC.se), "\n"))
  cat(paste("BIAS = ", formatC(object$BIAS, digits = 4), "\n"))
  cat(paste("Odds Ratio = ", formatC(object$theta, digits = 4), "\n"))
  cat(paste("Log Odds Ratio = ", formatC(object$log.theta, digits = 4), "\n"))
  cat(paste("Std. Err. for log Odds Ratio = ", formatC(object$LOR.se), "\n"))
    cat(paste("Odds Ratio Skill Score = ", formatC(object$orss, digits = 4), "\n")) 
  cat(paste("Std. Err. for Odds Ratio Skill Score = ", formatC(object$ORSS.se), "\n"))
    
  cat(paste("Extreme Dependency Score (EDS) = ", formatC(object$eds, digits = 4), "\n"))
  cat(paste("Std. Err. for EDS = ", formatC(object$eds.se, digits=4), "\n"))
  cat(paste("Symmetric Extreme Dependency Score (SEDS) = ", formatC(object$seds, digits = 4), "\n"))
  cat(paste("Std. Err. for SEDS = ", formatC(object$seds.se, digits=4), "\n"))
  cat(paste("Extremal Dependence Index (EDI) = ", formatC(object$EDI, digits=4), "\n"))
  cat(paste("Std. Err. for EDI = ", formatC(object$EDI.se, digits=4), "\n"))
  cat(paste("Symmetric Extremal Dependence Index (SEDI) = ", formatC(object$SEDI, digits=4), "\n"))
  cat(paste("Std. Err. for SEDI = ", formatC(object$SEDI.se, digits=4), "\n"))
  }

