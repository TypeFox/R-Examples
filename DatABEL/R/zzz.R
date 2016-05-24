.onAttach <- function(lib, pkg) {
  pkgVersion <- "0.9-6"
  welcomeMessage <- paste0(pkg, " v.", pkgVersion," loaded\n")

  .Call("checkNumBits");
  packageStartupMessage(welcomeMessage)
}
