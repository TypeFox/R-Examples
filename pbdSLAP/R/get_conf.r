### This file is only called by
###   "pbd*/src/Makevars.in" and "pbd*/src/Makevar.win"
### to find the default configurations from
###   "pbd*/etc${R_ARCH}/Makconf".
get.conf <- function(arg, arch = '', package = "pbdMPI"){
  file.name <- paste("./etc", arch, "/Makeconf", sep = "")
  file.path <- tools::file_path_as_absolute(
                 system.file(file.name, package = package))
  ret <- scan(file.path, what = character(), sep = "\n", quiet = TRUE)

  id <- grep(paste("^", arg, " = ", sep = ""), ret)
  if(length(id) > 0){
    cat(gsub(paste("^", arg, " = (.*)", sep = ""), "\\1", ret[id[1]]))
  } else{
    stop("The arg is not found.")
  }

  invisible()
} # End of get.conf().
