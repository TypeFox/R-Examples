
  cda <- new( "Module" )
  cd <- new( "Module" )
  array <- new( "Module" )
  dispersion <- new( "Module" )


 .onLoad <- function(libname, pkgname){
     loadRcppModules(direct=FALSE)
 }

