`.onLoad` <- 
function(libname, pkgname) {
}


`.onAttach` <- 
function(libname, pkgname) {
	if (interactive()) {
		packageStartupMessage('randomNames ',paste(paste(unlist(strsplit(as.character(packageVersion("randomNames")), "[.]")), c(".", "-", ".", ""), sep=""), collapse=""),'  For help type: help("randomNames")')
	}
}
