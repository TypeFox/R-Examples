# first and last
.onLoad <- function(libname,pkgname)
{
#require("utils",quietly=TRUE,warn.conflicts=FALSE)
#require("stats",quietly=TRUE,warn.conflicts=FALSE)
#require("gam",quietly=TRUE,warn.conflicts=FALSE)
#require("splines",quietly=TRUE,warn.conflicts=FALSE)
ver <- utils:::packageDescription(pkgname,fields="Version")
packageStartupMessage(paste("This is 'mtsdi' library",ver,sep=" "))
#library.dynam("mtsdi",pkg,lib)
}


.onUnload <- function(libpath)
{
#library.dynam.unload("mtsdi",libpath)
}

