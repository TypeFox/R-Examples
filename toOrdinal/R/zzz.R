`.onLoad` <- 
function(libname, pkgname) {
}


`.onAttach` <- 
function(libname, pkgname) {
	if (interactive()) {
		packageStartupMessage('toOrdinal ',paste(paste(unlist(strsplit(as.character(packageVersion("toOrdinal")), "[.]")), c(".", "-", ".", ""), sep=""), collapse=""),'  For help type: help("toOrdinal")')
	}
}
