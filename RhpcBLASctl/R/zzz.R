.onLoad <- function (lib, pkg) {
  library.dynam("RhpcBLASctl", pkg, lib)
}
