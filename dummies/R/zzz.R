
.onAttach <- function( libname, pkgname ) {

  packageStartupMessage( 
    pkgname ,
    "-" ,
    utils::packageVersion(pkgname, libname),
    " provided by Decision Patterns\n" ,
    domain = NA
  )
  
  # THIS SETS THE DEFAULTS FOR dummy.classes
  if( is.null( getOption("dummy.classes") ) )
    options( "dummy.classes" = c("factor","character") )

}


