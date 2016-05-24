.onLoad <- function(libname = find.package("mvQuad"), pkgname = "mvQuad") {
  utils::data("QuadRules", package=pkgname, envir=parent.env(environment()))
}
