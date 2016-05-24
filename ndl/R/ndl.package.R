#This is a general config section
ndl.package <- function() {
  library(help=ndl)$info[[1]] -> version
  if (!is.null(version)) {
    version <- version[pmatch("Version",version)]
    um <- strsplit(version," ")[[1]]
    version <- um[nchar(um)>0][2]
    hello <- paste("This is ndl version ",version,". \nFor an overview of the package, type 'help(\"ndl.package\")'.",sep="")
    packageStartupMessage(hello)
  } else {
    packageStartupMessage("Are we in devtools mode?")
  }
}

#set.ndl.options <- function()
## function used to set options
#{ 
  #add options here
  #ex: options()
#}

.onLoad <- function(...) {
#  set.ndl.options()
}

.onAttach <- function(...) { 
  ndl.package()
#  set.ndl.options()
}

.onUnload <- function(libpath) library.dynam.unload("ndl", libpath)
