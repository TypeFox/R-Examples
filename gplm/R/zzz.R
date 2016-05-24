.onLoad <- function(lib, pkg){
  library.dynam("gplm", pkg, lib)
  ## packageStartupMessage("gplm installed")
}
.onUnload <-function(libpath)
  library.dynam.unload("gplm", libpath)

