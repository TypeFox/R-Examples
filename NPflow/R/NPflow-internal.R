#'@keywords internal
.onUnload <- function (libpath) {
  library.dynam.unload("NPflow", libpath)
}

