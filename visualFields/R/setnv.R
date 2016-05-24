setnv <- function( nvtxt = "nvsapdefault" ) {

  nvinternal        <- eval( parse( text = nvtxt ) ) 
  nvinternal$nvname <- nvtxt
  assign( "nv", nvinternal, envir = visualFields::vfenv )

}
