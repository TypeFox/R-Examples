## new version
.onAttach <- function (libname, pkgname) {
  if (!interactive()) 
    return()
  Rcmdr <- options()$Rcmdr
  plugins <- Rcmdr$plugins
  if (!pkgname %in% plugins) {
    Rcmdr$plugins <- c(plugins, pkgname)
    options(Rcmdr=Rcmdr)
    if("package:Rcmdr" %in% search()) {
      if(!getRcmdr("autoRestart")) {
        options(Rcmdr=Rcmdr)
        closeCommander(ask=FALSE, ask.save=TRUE)
        Commander()
      }
    }
    else {
      Commander()
    }
  }
}

## older version
## .onAttach <- function(libname, pkgname){
##   if (!interactive()) return()
##   Rcmdr <- options()$Rcmdr
##   plugins <- Rcmdr$plugins
##   if ((!pkgname %in% plugins) && !getRcmdr("autoRestart")) {
##     Rcmdr$plugins <- c(plugins, pkgname)
##     options(Rcmdr=Rcmdr)
##     closeCommander(ask=FALSE, ask.save=TRUE)
##     Commander()
##   }
## }
