`%!in%` <- function(x,table) !`%in%`(x,table)

`%without%` <- function(x,table){
  x[ !x %in% table ]
}

`%x=%` <- function(txt, n){
	strrep( txt, n = n ) 
}

`%x=|%` <- function(txt, length.out){
  strrep( txt, length.out = length.out )
}

strrep <- function( txt, n, length.out = getOption("width") ){
	if( !missing( n ) ){
		paste( rep( txt, n ), collapse = "" )
	} else {
		out <- paste( rep(txt, ceiling( length.out / nchar(txt)  )) , collapse = "" )
		substr( out, 1, length.out)
	}
}

`%of%` <- function(x,y){
  inherits(x,y)
}

