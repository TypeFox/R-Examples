.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname);
  pkg <- AromaAffymetrix(pkgname);
  assign(pkgname, pkg, envir=ns);
} # .onLoad()


.onAttach <- function(libname, pkgname) {
  pkg <- get(pkgname, envir=getNamespace(pkgname));
  .setupAromaAffymetrix(pkg);
  startupMessage(pkg);
}
