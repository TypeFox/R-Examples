## covr: skip=all

.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname)
  pkg <- Package(pkgname)
  assign(pkgname, pkg, envir=ns)
} # .onLoad()


.onAttach <- function(libname, pkgname) {
  pkg <- get(pkgname, envir=getNamespace(pkgname))
  # require("affxparser") - install if missing
  .requireBiocPackage("affxparser", neededBy=pkgname)
  startupMessage(pkg)
}


############################################################################
# HISTORY:
# 2013-01-05
# o Added .onLoad().
# 2012-04-01
# o Now .onAttach() hides the fact that it tries to load 'affxparser'
#   from 'R CMD check'.
############################################################################
