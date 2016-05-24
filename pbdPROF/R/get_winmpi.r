### These functions are only called by "pbdMPI/src/Makevars.win" to obtain
### possible MPI dynamic/static library from the environment variable
### "MPI_ROOT" preset by users.

### Get conversion.
get.sysenv <- function(flag){
  sysenv <- Sys.getenv(flag)
  # sysenv <- shortPathName(sysenv)
  sysenv <- gsub("\\\\", "/", sysenv)
  if(sysenv != ""){
    if(length(grep("/$", sysenv)) == 0){
      sysenv <- paste(sysenv, "/", sep = "")
    }
  }
  cat(sysenv)

  invisible()
} # End of get.sysenv().

