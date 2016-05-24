##################################################
###                Global.R
setGeneric("privateA",function(object){standardGeneric("privateA")})
setGeneric("publicA",function(object){standardGeneric("publicA")})
setGeneric("publicB",function(objectV,objectW){standardGeneric("publicB")})
functionClassicA <- function(age){return(age*2)}

#.onLoad <- function(libname,pkgname){packageStartupMessage("### [packS4] You load package <",pkgname,"> from <",libname,"> ###\n",sep="")}
#.onAttach <- function(libname,pkgname){packageStartupMessage("### [packS4] You attach pacakge <",pkgname,"> from <",libname,"> ###\n",sep="")}
.onUnload <- function(libpath){cat("### [packS4] You unload the pacakge ###\n",sep="")}

