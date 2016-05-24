.onAttach <- function( libname, pkgname ) {

  packageStartupMessage( 
    pkgname ,
    "-" ,
    utils::packageVersion(pkgname, libname),
    " - Copyright \u00a9 ", substr(Sys.Date(),1,4),
    " Decision Patterns" ,
    domain = NA
  )

}

