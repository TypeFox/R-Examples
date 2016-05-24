##' R modules for JASPAR databases
##' 
##'
##' \tabular{ll}{
##' Package: \tab JASPAR \cr
##' Type: \tab Package \cr
##' License: \tab GPL (>= 2) \cr
##' }
##' 
##' This package contains: \cr
##' modules and functions for JASPAR data processing and visualization.
##' \cr
##' 
##' An overview of functions
##' \tabular{ll}{
##' Function \tab Description\cr
##' \code{make_template} \tab  Make a template that feeds into JASPAR databases
##' }
##' 
##'
##' ## To install from online repositories (e.g. CRAN)  \cr 
##' install.packages(pkgs="JASPAR", repos="http://cran.r-project.org") \cr
##' 
##' ## Sometimes the offical repository might not be up to date, \cr
##' ## then you may install it from a downloaded source file; \cr
##' ## please replace '<current-version>' with actual version numbers: \cr
##' install.packages(pkgs="JASPAR_<current-version>.tar.gz", repos=NULL) \cr
##' 
##' 
##' ## Load the package and get a complete list of functions, use \cr 
##' library(JASPAR) \cr 
##' ls("package:JASPAR") \cr
##'
##' ## help documantation of the package \cr 
##' help(JASPAR)  # this page
##' 
##' @name JASPAR-package
##' @aliases JASPAR JASPAR-package
##' @docType package
##' @title The Package for R modules for JASPAR databases
##' @references Zhao et al (2013), 
##' "JASPAR 2013: An extensively expanded and updated open-access database of
##' transcription factor binding profiles." (\emph{In preparation})
##' 
##'
##'
##' See \code{citation("JASPAR")} for BibTeX entries for LaTeX users.
##' @author Xiaobei Zhao
##' 
##' Maintainer: Xiaobei Zhao \email{xiaobei (at) binf.ku.dk}
##' @concept JASPAR 
##' @keywords package
##' @return \code{NULL}
##' @seealso \code{\link{make_template}}
##' @examples
##' require(JASPAR)        # load JASPAR
##' help(JASPAR)           # a help document of JASPAR 
##' ## data(package="JASPAR") # a list of datasets available in JASPAR (TBA)
##' ls("package:JASPAR")   # a list of functions available in JASPAR
##' help(package="JASPAR") # help documentation on JASPAR
##' citation("JASPAR")     # citation for publications
##' demo("JASPAR-demo")    # run the demo 
##'
##' ## view JASPAR Description
##' packageDescription("JASPAR")
##' 
##' ## ## view JASPAR vignette (TBA)
##' ## vignette("JASPAR-vignette",package="JASPAR")
##' 
NULL
