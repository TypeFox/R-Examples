### send to file

..redirect <- function(object,file,type = c("output", "message", "both"),append){
	type <- match.arg( type )
  if( is.character(file) ) file <- file( file, open = if(append) "at" else "wt")
  if( type %in% c("output", "message") ) sink( file , append = append, type = type)
  else{
     sink( file, type = "output")
     sink( file, type = "message")
  }
	# borrow the capture function in svMisc instead of all that
	# or should I create a text connection, write the expression in, and source it ?
  txt <- capture.output( 
	  eval( parse( text = deparse(substitute(object))) )
		)
	if( type %in% c("output", "both") ) cat( txt, sep = "\n" )	
	if( type %in% c("output", "message") ) sink(type = type)
  else{            
     sink( type = "output")
     sink( type = "message")
  }
  close(file)
	
}

..readfromfile <- function(object, file, append = FALSE, objname, envir, verbose = getOption("verbose")){
  content <- readLines(file,n=-1)
  if(append){ # try to get the content of the object
    obj.content <- try( get(objname,envir), silent = !verbose )
    if( class(obj.content) == "try-error" ) obj.content <- NULL
    content <- c(obj.content,content)
  }
  assign(objname,content, envir=envir )   
  invisible(NULL)
}

### send to file
`%>%`   <- function( object, file) ..redirect(object,file,type="output" , append=FALSE)
`%>>%`  <- function( object, file) ..redirect(object,file,type="output" , append=TRUE )
`%2>>%` <- function( object, file) ..redirect(object,file,type="message", append=TRUE )
`%2>%`  <- function( object, file) ..redirect(object,file,type="message", append=FALSE)
`%*>>%` <- function( object, file) ..redirect(object,file,type="both", append=TRUE )  
`%*>%`  <- function( object, file) ..redirect(object,file,type="both", append=FALSE)    
  
### read from file
`%<%`   <- function( object, file) ..readfromfile(object,file, append=FALSE, objname=deparse(substitute(object)), envir = parent.frame(1))
`%<<%`  <- function( object, file) ..readfromfile(object,file, append=TRUE , objname=deparse(substitute(object)), envir = parent.frame(1))

`%|%` <- function( r, u){
  tf <- tempfile()
  on.exit( unlink( tf) )
  sink(tf)
  print( r )
  sink()
  pi <- pipe( paste( "cat", tf, "|", u ) )
  out <- readLines( pi )
  close(pi)
  structure( out, class= "unixoutput")
}

print.unixoutput <- function( x, ...){
  if( length(x) ){
	cat( x, sep = "\n", ...)
  }
}

