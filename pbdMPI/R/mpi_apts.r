### Obtain MPI array pointers.

arrange.mpi.apts <- function(){
  rm(list = c(".__MPI_APTS__"), envir = .GlobalEnv)
  if(is.loaded("spmd_initialize", PACKAGE = "pbdMPI")){
    .Call("arrange_MPI_APTS", PACKAGE = "pbdMPI")
  } else{
    stop("The MPI daemon may be stop.")
  }
  invisible()
} # End of arrange.mpi.apts().

