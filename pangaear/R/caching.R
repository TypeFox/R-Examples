env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  path <- rappdirs::user_cache_dir("pangaear")
  env$path <- path
}
