##' .onAttach function
##'
##' A function to print a welcome message on loading package  
##'
##' @param libname libname argument
##' @param pkgname pkgname argument
##' @return ...
##' @export

.onAttach <- function(libname, pkgname)
{
	packageStartupMessage("\n Welcome to 'spatsurv': Spatial Survival Analysis\n B. M. Taylor & B. S. Rowlingson.\n
Type 'spatsurvVignette()' to view the package vignette.\n
Type 'citation(\"spatsurv\")' to view the citation for this package.\n
Please see the spatsurv package NEWS file for latest additions, changes and bug fixes.", appendLF=T)
}

