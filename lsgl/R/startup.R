.onAttach <- function(libname, pkgname) {
	
	c.config <- lsgl.c.config()
	
	if(c.config$debugging) packageStartupMessage("lsgl: Compiled with debugging on -- this may slow down runtime")
	
}