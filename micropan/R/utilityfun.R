#$Id: utilityfun.R 142 2013-07-23 12:04:48Z larssn $

gregexpr <- function( pattern, text, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, extract=FALSE ){
  lst <- base::gregexpr( pattern, text, ignore.case, perl, fixed, useBytes )
  if( extract ){
    lst <- lapply( 1:length( lst ), function( i ){ substring( text[i], lst[[i]], lst[[i]]+attr( lst[[i]], "match.length" ) - 1 ) } )
  }
  return( lst )
}



# isIn <- function( patterns, text, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE, invert = FALSE, cl = NULL  ){
#   npat <- length( patterns )
#   if( is.null( cl ) ){
#     isin <- unlist(lapply(lapply(1:npat, function(i) grep( patterns[i], text, ignore.case, perl, fixed, useBytes, invert )),length))>0
#   } else {
#     text <- eval(text)
#     isin <- unlist(lapply(parLapply(cl,1:npat, function(i) grep( patterns[i], text, ignore.case, perl, fixed, useBytes, invert )),length))>0
#   }
#   return( isin )
# }


