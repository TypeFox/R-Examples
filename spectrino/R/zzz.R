.spnEnv <- new.env()

.onLoad <- function(libname, pkgname) {
  assign(".sprnOpt", list("timeOut", "host", "port"), envir = .spnEnv) 
  assign(".wsConnected", FALSE, envir = .spnEnv)  
  assign(".wsSprnFlag", FALSE, envir = .spnEnv) 
  assign(".wsSprnText", "", envir = .spnEnv)
  assign(".wsSprnHandle", 0, envir = .spnEnv)
  assign(".wsSprn", 0, envir = .spnEnv)
}