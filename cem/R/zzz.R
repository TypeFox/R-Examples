.noGenerics <- TRUE

.onAttach <- function(libname, pkgname){
  packageStartupMessage("\nHow to use CEM? Type vignette(\"cem\")\n")
}
