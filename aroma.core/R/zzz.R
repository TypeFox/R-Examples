.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname);
  pkg <- AromaCore(pkgname);
  assign(pkgname, pkg, envir=ns);

  .loadSettings(pkgname)
} # .onLoad()


.onAttach <- function(libname, pkgname) {
  pkg <- get(pkgname, envir=getNamespace(pkgname));
  .setupAromaCore(pkg);
  startupMessage(pkg);
}
