# replace some spaces in a sentence by newlines to get a better looking
# cloud


#' Replace some spaces in multi-word sentences by newlines
#' 
#' Replace a space character by a newline in a multi-word sentence to get a
#' better height / width ratio
#' 
#' Very long tags, for example GO Term descriptions, make a bad tag cloud.
#' \code{strmultline} tries to chop up such a long sentence into multiple
#' (currently two) lines, to get a better height / width ratio.
#' 
#' @param strings a character vector containing the multi-word sentences to be
#' split
#' @param ratio the desired ratio height / width
#' @return A character vector containing the modified sentences.
#' @author January Weiner <january.weiner@@gmail.com>
#' @seealso \code{\link{tagcloud}}
#' @keywords strings splitting
#' @export strmultline
strmultline <- function( strings, ratio= 0.2 ) {

  strings <- as.character( strings )
  n <- length( strings )
  splits <- strsplit( strings, "[[:space:]]" ) 

  for( i in 1:n ) {

    x <- splits[[i]]

    nw <- length( x )
    if( nw == 1 ) next

    sl <- sapply( x, nchar )
    r1 <- 1 / nchar( strings[i] )

    if( r1 > ratio ) next 

    r2s <- c()

    # try all possible splits. Surely there is a more effective approach,
    # feel free to modify the code.
    for( j in 1:(nw - 1) ) {

      str.n1 <- paste( x[1:j], collapse= " " )
      str.n2 <- paste( x[(j+1):nw], collapse= " " ) 

      r2 <- 2 / max( nchar( str.n1 ), nchar( str.n2 )  )

      r2s <- c( r2s, r2 )
    }

    j <- which.min( abs( r2s - ratio ) )

    str.n <- paste( paste( x[1:j], collapse= " " ), paste( x[(j+1):nw], collapse= " " ), sep= "\n" )
    strings[i] <- str.n 
  } 
 
  return( strings )
}
