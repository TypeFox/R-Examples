### These functions are only called by "pbdNCDF4/src/Makevars.win" to obtain
### possible NCDF4 static library from the environment variable
### "NCDF4_ROOT_64", "NCDF4_ROOT_32", and "NCDF4_ROOT" preset by users.

### Get conversion.
get.sysenv <- function(flag){
  sysenv <- Sys.getenv(flag)
  sysenv <- shortPathName(sysenv)
  sysenv <- gsub("\\\\", "/", sysenv)
  if(sysenv != ""){
    if(length(grep("/$", sysenv)) == 0){
      sysenv <- paste(sysenv, "/", sep = "")
    }
  }
  cat(sysenv)

  invisible()
} # End of get.sysenv().

