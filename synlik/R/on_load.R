## preliminary code

.printSynlikVersion <- function()
{ library(help=synlik)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  cat(paste("This is \"synlik\" ",version,"\n"),sep="")
}


.onAttach <- function(...) { 
  .printSynlikVersion()
}

.onUnload <- function(libpath) library.dynam.unload("synlik", libpath)

.onLoad <- function(lib,pkg) {
   library.dynam("synlik", pkg, lib)
}
