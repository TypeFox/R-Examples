`.onLoad` <- function(libname, pkgname) {
  assign( "sd.default", stats::sd , asNamespace(pkgname))
  assign("var.default", stats::var, asNamespace(pkgname))
}
