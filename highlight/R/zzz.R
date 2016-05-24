#' @importFrom grDevices rgb postscript dev.off pdf
#' @importFrom tools vignetteEngine
#' @useDynLib highlight
NULL

NAMESPACE <- environment()

.onLoad <- function(libname, pkgname){
	options( detective = simple_detective )
	
	vignetteEngine("highlight", weave = Hweave, tangle = Htangle,
        pattern = "[.][hHrRsS]nw$", 
        package = "highlight")
}

