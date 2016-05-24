`.onLoad` <- 
function(libname, pkgname) {
}


`.onAttach` <- 
function(libname, pkgname) {
	if (interactive()) {
		packageStartupMessage('SGPdata ',paste(paste(unlist(strsplit(as.character(packageVersion("SGPdata")), "[.]")), c(".", "-", ".", ""), sep=""), collapse=""),' (3-1-2016).  For help: >help("SGPdata")')
	}
}
