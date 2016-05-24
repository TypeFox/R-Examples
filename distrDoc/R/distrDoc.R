.onLoad <- function(lib, pkg) { # extended 03-28-06: P.R.
#    require("methods", character = TRUE, quietly = TRUE)
}

.onAttach <- function(libname, pkgname)
{                                                                    
#    if (.Platform$OS.type == "windows" && require("Biobase")
#        && interactive() && .Platform$GUI == "Rgui")
#        addVigs2WinMenu("distrDoc")

buildStartupMessage(pkg="distrDoc",  library=libname, 
                    packageHelp=TRUE, # MANUAL="http://www.uni-bayreuth.de/departments/math/org/mathe7/DISTR/distr.pdf",
                    VIGNETTE=gettext("This package provides a vignette to package \"distr\" and to several related packages; try vignette(\"distr\")."))
  invisible()
} 

### just to have it -- not to export it
#setClass("Integer", contains ="numeric",
#          validity = function(object) all(.isInteger(object)))
