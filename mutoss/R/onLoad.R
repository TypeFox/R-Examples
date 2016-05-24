.onLoad <- function(libname, pkgname) {
#	if(!require("multtest", character.only=TRUE)) {
#		if (interactive()) {
#			answer <- readline("Multtest is missing - do you want to install it (y/N)? ")
#			if (substr(answer, 1, 1) %in% c("y","Y")) {
#				source("http://bioconductor.org/biocLite.R")
#				biocLite("multtest")
#				require("multtest")
#			} else {
#				warning("Package multtest is not avaible - please install it!")
#			}
#		} else {
#			warning("Package multtest is not avaible - please install it!")
#		}
#	}
}  
