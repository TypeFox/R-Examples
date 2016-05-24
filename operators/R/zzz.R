.onAttach <- function( libname, pkgname ){
  loadOptions( )
}

loadOptions <- function( ) {
	opt.file <- system.file( "options", "options.R", package = "operators" )
	if( file.exists( opt.file ) ){
		source( opt.file )
	}
}

