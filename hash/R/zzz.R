
# .First.lib <- function( libname, pkgname ) {
.onAttach <- function( libname, pkgname ) {

  packageStartupMessage( 
    pkgname ,
    "-" ,
    utils::packageVersion(pkgname, libname),
    " provided by Decision Patterns\n" ,
    domain = NA
  )

}

