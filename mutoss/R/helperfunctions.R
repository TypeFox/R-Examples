
requireLibrary <- function(package) {
	if(!require(package, character.only=TRUE)) {
		answer <- readline(paste("Package ",package," is required - should we install it?",sep=""))
		if (substr(answer, 1, 1) %in% c("y","Y")) {
			install.packages(package)
			require(package, character.only=TRUE)
		} else {
			stop(paste("Required package",package,"should not be installed"))
		}
	}
}

