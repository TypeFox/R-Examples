# (c) William R. Engels, 2015

#' @importFrom utils packageDescription
#' @importFrom stats dchisq pchisq qchisq runif

.onAttach <- function(libname, pkgname){
	version <- packageDescription("HWxtest", fields = "Version")
	packageStartupMessage("------------------------------------------------")
	packageStartupMessage("                  HWxtest")
	packageStartupMessage(paste("               version", version))
	packageStartupMessage("Please see the package vignette for instructions")
	packageStartupMessage("------------------------------------------------")
}