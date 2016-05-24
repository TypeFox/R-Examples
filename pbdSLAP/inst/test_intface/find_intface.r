### This file is modified from "BRugs/R/windows/findOpenBUGS.R".

### This function is only called by "RmpiSPMD/src/Makevars.in" to obtain
### possible BLACS INTFACE for Fortran to C..
find.intface <- function(){
  system("./00_make.sh", intern = TRUE, ignore.stderr = TRUE)
  invisible()
} # End of find.intface().

