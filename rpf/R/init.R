.onLoad <- function(libname, pkgname) {
	if ("package:parallel" %in% search()) {
		cores <- detectCores(logical=TRUE)
		.Call(setNumberOfCores, cores)
		#packageStartupMessage(paste("OpenMP will use", cores, "cores"))
	}
}

.onAttach <- function(libname, pkgname) {
	if (! .Call(hasOpenMP_wrapper)) {
		packageStartupMessage("RPF is not compiled to take advantage of computers with multiple cores.")
	}
}
