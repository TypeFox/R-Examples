# Ensure that the DLL is unloaded when unloading this library.
# Get a list of loaded DLLs with getLoadedDLLs()
.onUnload <- function(libpath) {
  library.dynam.unload("autovarCore", libpath)
}
