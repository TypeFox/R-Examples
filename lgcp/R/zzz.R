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
	packageStartupMessage("\n Welcome to 'lgcp': Log-Gaussian Cox Process\n B. M. Taylor & T. M. Davies & B. S. Rowlingson & P. J. Diggle.\n
Type '?lgcp' for details.\n
Type 'lgcpvignette()' to view the basic package vignette.\n
Type 'lgcpbayes()' to view the Bayesian package vignette.\n
Type 'citation(\"lgcp\")' to view the citation for this package.\n
Please see the lgcp package NEWS file for latest additions, changes and bug fixes.)", appendLF=T)
}

