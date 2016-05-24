### Lastest load into a package.

### Export Namespace does not use .First.lib() and .Last.lib(), but use
### .onLoad() and .onUnload().

.onLoad <- function(libname, pkgname){
  if(! is.loaded("spmd_initialize")){
    library.dynam("pbdMPI", "pbdMPI", libname)
    if(pbdMPI::comm.is.null(0L) == -1){
      pbdMPI::init()
    }
  }

  if(! is.loaded("slap_blacs_gridinit")){
    library.dynam("pbdSLAP", "pbdSLAP", libname)
  }

  if(! is.loaded("R_blacs_init")){
    library.dynam("pbdBASE", "pbdBASE", libname)
  }
  
  pbdMPI::pbd_opt(BLDIM=c(16, 16))
  pbdMPI::pbd_opt(ICTXT=0)

  invisible()
}

.onUnload <- function(libpath){
  library.dynam.unload("pbdDMAT", libpath)
  invisible()
}
