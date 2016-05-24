### This file is only called by
###   "pbd*/src/Makevars.in" and "pbd*/src/Makevar.win"
### to find the default configurations from
###   "pbd*/etc${R_ARCH}/Makconf".
get.lib <- function(arg, arch, package = "pbdPROF"){
  file.name <- paste("./etc", arch, "/Makeconf", sep = "")
  file.path <- tools::file_path_as_absolute(
                 system.file(file.name, package = package))
  ret <- scan(file.path, what = character(), sep = "\n", quiet = TRUE)

  id <- grep(paste("^", arg, " = ", sep = ""), ret)
  if(length(id) > 0){
    lib.file <- gsub(paste("^", arg, " = (.*)", sep = ""), "\\1", ret[id[1]])
    file.path <- tools::file_path_as_absolute(
                   system.file(lib.file, package = package))
    cat(file.path)
  } else{
    comm.stop("The arg is not found.")
  }

  invisible()
} # End of get.lib().
