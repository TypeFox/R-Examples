### Lastest load into a package.

### Export Namespace does not use .First.lib() and .Last.lib(), but use
### .onLoad() and .onUnload().
# .First.lib <- function(lib, pkg){
# } # End of .First.lib().

# .Last.lib <- function(libpath){
# } # End of .Last.lib().

.onLoad <- function(libname, pkgname){
  if(! is.loaded("spmd_initialize", PACKAGE = "pbdMPI")){
    library.dynam("pbdMPI", "pbdMPI", libname)
    if(pbdMPI::comm.is.null(0L) == -1){
      pbdMPI::init()
    }
  }

  library.dynam("pbdSLAP", pkgname, libname)
  invisible()
} # End of .onLoad().

.onUnload <- function(libpath){
  pbdSLAP::slap.finalize()
  library.dynam.unload("pbdSLAP", libpath)
  invisible()
} # End of .onUnload().

