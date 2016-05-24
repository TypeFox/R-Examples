#' On attach
#'
#' @author Tomislav Hengl \email{tom.hengl@wur.nl} and Andrew Sila \email{asila@cgiar.org}

.onLoad <- function(libname, pkgname) {
  data("wavenumbers", package=pkgname, envir=parent.env(environment()))
}

.onAttach <- function(libname, pkgname)  {
  
  ## print on start-up:
	# pkg.info <- utils::packageDescription('GSIF')
	pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION", package=pkgname, lib.loc=libname), fields=c("Version","Date")))
	packageStartupMessage(paste(pkgname, " version ", pkg.info["Version"], " (", pkg.info["Date"], ")", sep=""))

  spec.env()
  packageStartupMessage(paste("URL: http://worldagroforestry.org"))

  return(invisible(0))
  	
}

# end of script;