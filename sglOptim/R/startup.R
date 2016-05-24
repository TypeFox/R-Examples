.onAttach <- function(libname, pkgname) {
	
	c.config <- sgl.c.config()
	
	if(!c.config$omp.supported) packageStartupMessage("sglOptim: openMP (multithreading) is not supported on this system")
	if(c.config$debugging) packageStartupMessage("sglOptim: Compiled with debugging on -- this may slow down runtime")
	
}