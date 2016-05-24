.onLoad <- function(libname, pkgname) {
  # at the compile time a flag is set which can
  # be accessed by using the function isDllVersion 
  if ( !.Call( "isDllVersion") ){
      pyConnect()
  } else if ( nchar(Sys.getenv('PYTHON_EXE')) > 0 ) {
  	  pyConnect(pythonExePath=Sys.getenv('PYTHON_EXE'))
  }
  
  invisible(NULL)
}

