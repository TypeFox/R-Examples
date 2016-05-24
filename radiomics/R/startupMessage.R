#' @useDynLib radiomics
#' @importFrom Rcpp sourceCpp

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("radiomics: Functions for texture analysis of greyscale images.\nEnter ?calc_features to see available texture features.")
}

.onUnload <- function (libpath) {
  library.dynam.unload("radiomics", libpath)
}