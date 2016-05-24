.onLoad <- function(libname, pkgname) {

 # DEFINE AUTOHELP ACTIVE BINDING
   # makeActiveBinding( "opt_help" , opt_help, baseenv() )
   
   # set optigrab options
   options( optigrab = list( 
       help = list()
     , style = gnu_style
     , options = list() 
     )
   )

}


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

