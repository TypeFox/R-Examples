setClass("rose",
         representation(rho="matrix",
                        cyclVar="numeric",
                        circle="numeric"),
         validity=val.rose)

setMethod("plot", signature(x = "rose",y = "missing"), plot.rose)

.onLoad <- function(libname, pkgname){
  library.dynam("IDPmisc", package=pkgname, lib.loc=libname)
}

.onUnLoad <- function(libpath){
  library.dynam.unload("IDPmisc",libpath=libpath)
}
