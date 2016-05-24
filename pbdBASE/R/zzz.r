### Lastest load into a package.

## Export Namespace does not use .First.lib() and .Last.lib(), but use
## .onLoad() and .onUnload().

#' @export
.Last.lib <- function(libpath){
  ### To free all BLACS points.
  base.finalize(mpi.finalize = FALSE)
} # End of .Last.lib().

.onLoad <- function(libname, pkgname){
  if(! is.loaded("spmd_initialize", PACKAGE = "pbdMPI")){
    library.dynam("pbdMPI", "pbdMPI", libname)
    if(pbdMPI::comm.is.null(0L) == -1){
      pbdMPI::init()
    }
  }

  library.dynam("pbdBASE", pkgname, libname)
  invisible()
} # End of .onLoad().

.onUnload <- function(libpath){
  library.dynam.unload("pbdBASE", libpath)
  invisible()
} # End of .onUnload().

