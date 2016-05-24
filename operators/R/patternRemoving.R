
### pattern removing
`%-~%` <- function(txt, pattern){
  .gsub <- ..gsub %but% getOption("operators.gsub")
  .gsub( pattern , "", txt)
}

# filters and remove
`%-~|%` <- function(txt, pattern){
  (txt %~|% pattern) %-~% pattern
}

`%o~|%` <- function(txt, pattern){
  txt <- txt %~|% pattern
	txt %o~% pattern
}

`%o~%` <- function(txt, pattern){
	if( txt %!~+% pattern) return(NULL)
	if( pattern %!~% "\\(.*?\\)" ) {
		pattern <- sprintf("(%s)", pattern) 
	} 
	if( pattern %!~% "^\\^" ){
		pattern <- sprintf( "^.*?%s", pattern ) 
	}
	if( pattern %!~% "\\$$" ){
	  pattern <- sprintf( "%s.*?$", pattern)
	}
	
	# how many chunks to keep
	n <- length( gregexpr("\\([^)]*\\)", pattern)[[1]]  ) 
	
	out <- rep( list(NULL), n )
	for( i in 1:n ){
		out[[i]] <- ifelse( txt %~% pattern, 
			gsub( pattern, sprintf("\\%d", i), txt, perl = TRUE ), 
			getOption("operators.o.nomatch") )
	}
	out <- do.call( cbind, out )
	rownames( out ) <- txt
	out
}



`%/~%` <- function( txt, rx ){
  .strsplit <- strsplit %but% getOption("operators.strsplit")
  unlist( .strsplit( txt, rx) )
}

`%s~%` <- function( txt, pattern ){
  if( pattern %!~% "^/") stop( gettext("the regular expression should start with a '/'") )
  pattern <- ( pattern %/~% "/" ) [-1]
  modif <- if( length(pattern) ==3 && nchar(pattern[3]) > 0 ){ # get the modifiers
    pattern[3]
  } else  getOption("operators.gsub")
  
  .gsub <- ..gsub %but% modif
  .gsub( pattern[1], pattern[2], txt )
}

### gsub or sub depending on the global argument
..gsub <- function(pattern, replacement, x, ignore.case = FALSE, 
    perl = FALSE, fixed = FALSE, useBytes = FALSE, global=TRUE){
  
   if(global) gsub(pattern = pattern ,replacement = replacement, x = x, 
   	ignore.case = ignore.case, perl = perl , fixed = fixed, useBytes = useBytes)
   else sub(pattern = pattern,replacement = replacement, x = x, 
   	ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
}


