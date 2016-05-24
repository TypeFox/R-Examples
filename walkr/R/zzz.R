
## CMD Check requires that the compressed vignette from building
## the package be small in size. Thus, we must compress it 
## in the zzz function

.onLoad <- function(libname, pkgname) {
  
  ## compress the vignette PDF to fix CMD Check WARNING
  
  Sys.setenv("_R_BUILD_COMPACT_VIGNETTES_" = "gs+qpdf")
}