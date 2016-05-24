#' dynamic loading
#' @param lib
#' @param pkg
#' @import RcppProgress
#' @author	karl

.onLoad <- function(lib, pkg) {
	library.dynam("RcppProgress", pkg, lib )
}

.onUnload <- function (libpath) {
	library.dynam.unload("RcppProgress", libpath)
}