.onAttach <- function(libname, pkgname) {
    actual <- utils::packageDescription(pkgname)[["Version"]]

    packageStartupMessage(paste("Welcome in package", pkgname, "version", actual))
    
#	essai <- available.packages(contriburl = "http://cran.at.r-project.org/")
	
#	if (dim(essai)[1]==0) {
#	        m <- paste("Your version is ", actual)
#            m <- paste(m, ". No internet connection is available to check for update.")
#            packageStartupMessage(m)

#	} else {
#		recent <- essai[pkgname, "Version"]
#	    if (as.numeric(actual) < as.numeric(recent)) {
#            m <- paste("Your version is ", actual, ". Most recent is ", recent, sep="")
#		    packageStartupMessage(m)
#            packageStartupMessage("Use update.packages() to update...")
#        }

#	}

}
