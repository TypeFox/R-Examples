
.onLoad <- function(libname,pkgname)
{
ver <- packageDescription("pgam",fields="Version")
require("stats",character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE)
require("utils",character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE)
cat("This is pgam library version",ver,"\n",sep=" ")
library.dynam("pgam",pkgname,libname)
}


.onUnload <- function(libpath)
{
# unloading process 
library.dynam.unload("pgam",libpath)
}
