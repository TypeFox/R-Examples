.requirePkg <- function(name, quietly=FALSE) {
  # Nothing to do?
  if (is.element(sprintf("package:%s", name), search())) {
    return(invisible(TRUE));
  }
  if (quietly) {
    # Load the package (super quietly)
    res <- suppressPackageStartupMessages(require(name, character.only=TRUE, quietly=TRUE));
  } else {
    res <- require(name, character.only=TRUE);
  }
  if (!res) throw("Package not loaded: ", name);
  invisible(res);
} # .requirePkg()


.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname);
  pkg <- Package(pkgname);
  assign(pkgname, pkg, envir=ns);
} # .onLoad()


.onAttach <- function(libname, pkgname) {
  pkg <- get(pkgname, envir=getNamespace(pkgname));
  .setupAromaCn(pkg);
  startupMessage(pkg);
}
