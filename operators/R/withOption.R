
options <- function( ... ){
	structure( base::options( ... ), class = "options" )
}

print.options <- function( x, ... ){
	print( unclass( x ), ... )
}

with.options <- function( data, expr, ...){
	old.op <- data; on.exit( base::options( old.op ) )
	out <- eval( substitute( expr ) )
	attr( out, "withOptions" ) <- options()[names(old.op)]
	class( out ) <- c( "withOptions", class(out) )
	out
}

print.withOptions <- function( x, verbose = getOption("verbose"), ...){
	old.op <- do.call( base::options, attr(x, "withOptions" ) )
	on.exit(base::options(old.op))
	if(verbose) cat ( "with options", paste( "\n", names(old.op), ":", old.op) ) 
	class( x ) <- setdiff( class(x), "withOptions" )
	print( x, ... )
}

