
.onLoad <- function(lib,pkg){
	library.dynam("RFgroove", package=pkg, lib.loc=lib)
	return(invisible(0))
}
