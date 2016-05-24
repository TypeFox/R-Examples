#  zzz.R
#
# System functions for the AGD library

#------------------------------.onAttach-------------------------------
.onAttach <- function(...){
  d <- packageDescription("AGD")
  packageStartupMessage(paste(d$Package,d$Version,d$Date))
  return()
}

version <- function(pkg="AGD"){
  lib <- dirname(system.file(package = pkg))
  d <- packageDescription(pkg)
  return(paste(d$Package,d$Version,d$Date,lib))
}
