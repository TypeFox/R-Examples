## Copyright (C) 2013-2015 Kasper Kristensen
## License: GPL-2

## .First.lib <- function(lib, pkg) {
##   library.dynam("TMB", pkg, lib)
## }

.onLoad <- function(lib, pkg) {
  library.dynam("TMB", pkg, lib)
}

## .LastLib <- function(libpath)
## {
##   library.dynam.unload("TMB", libpath)
## }


.onAttach <- function(lib, pkg) {
  exfolder <- system.file("examples", package = "TMB")
  dll <- paste0(exfolder, Sys.getenv("R_ARCH"), "/simple", .Platform$dynlib.ext)
  if(!file.exists(dll)) runExample("simple", dontrun=TRUE)
}
