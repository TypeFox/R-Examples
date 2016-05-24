.onLoad <- function(libname, pkgname){
  #Set MAgClass verbosity level to 1 (warnings, but no notes)
  options(magclass.verbosity=1)
  #Due to a function name conflict with the grid library 
  #(both packages contain the function getNames), 
  #it has to be assured that the grid library is loaded before the magclass library
  #if the grid library is installed
  suppressWarnings(try(do.call(what="library",args=list("grid"))))  
}