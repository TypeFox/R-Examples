.onLoad <- function(lib, pkg) {
  library.dynam("isa2", pkg, lib, local=FALSE);

  ## Default ISA options
  isa.options[["verbose"]] <-  FALSE
  isa.options[["status.function"]] <- function(...) isa.status.default(...)
  
}

.onUnload <- function(libpath) {
  library.dynam.unload("isa2", libpath)
}

