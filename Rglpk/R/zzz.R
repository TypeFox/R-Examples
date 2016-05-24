.onLoad <- function(libname, pkgname){
  ## initialize print routines
  .C("Rglpk_initialize", PACKAGE = pkgname)
}

.onAttach <- function(libname, pkgname){
    version <- .C( "Rglpk_get_engine_version",
                   GLPK_version = character(1L), PACKAGE = pkgname )
    packageStartupMessage( sprintf("Using the GLPK callable library version %s", version) )
}
