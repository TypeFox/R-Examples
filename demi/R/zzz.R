#===========================================================================#
# xps package initialization
#===========================================================================#

.onAttach <- function(libname, pkgname) {
#	require(methods);
#	require(utils);
	msg <- paste("\nWelcome to", pkgname, "version", utils::packageDescription(pkgname, fields="Version"), "\n");
	msg <- paste( msg, "    An R wrapper for DEMI - estimating Differential Expression from Multiple Indicators\n", sep = "" )
	msg <- paste( msg, "    Easy to use tool for analysing differential expression on high-density microarray platforms\n", sep = "" )
	msg <- paste( msg, "    (c) Copyright 2013 by Sten Ilmjarv and Hendrik Luuk\n", sep = "" )
	msg <- paste( msg, "    \n", sep = "" )
#	cat(paste("\nWelcome to", pkgname, "version", packageDescription(pkgname, fields="Version"), "\n"));
#	cat("    DEMI - Differential Expression from Multiple Indicators\n");
#	cat("    An R package for analysing differential expression on high-density microarray platforms\n");
#	cat("    (c) Copyright 2013 by Sten Ilmjarv and Hendrik Luuk\n");
#	cat("    \n");
	packageStartupMessage(msg);
}

#.onUnload <- function(libpath) {
#	library.dynam.unload( "demi", libpath )
#}
