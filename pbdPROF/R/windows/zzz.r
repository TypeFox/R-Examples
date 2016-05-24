.onAttache <- function(libname, pkgname){
  msg <-
"
    No MPI profiler support for Windows system.
"
  packageStartupMessage(msg)
  invisible()
} # End of .onAttach().
