## ## .First.lib is meaningful for S-Plus only, not for R.
## ## It doesn't get defined in S-Plus when it is inside the if.R() statement
## ## .First.lib <- function(library, section) { ## S-Plus notation
## .First.lib <- function(libname, pkgname) {  ## R notation, to be used by S-Plus
##   invisible(options(HH.ROOT.DIR=paste(libname, pkgname, sep="/")))
##   attach(hh.new("HH.Data"))
## }

## if.R(r={
  .onLoad <- function(libname, pkgname) {
    options(HH.ROOT.DIR=paste(libname, pkgname, sep="/"))
   }
##   rm(.First.lib) ## not meaningful in R, hence removed
## }, s={}
##      )
