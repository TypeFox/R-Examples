.onLoad <- function( libname, pkgname ) {
    version <- detect_pandoc_version()$num
    if ( is.numeric(version) ) set_pandoc_version(version)
}