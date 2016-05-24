# minimalistic "editor" to allow interactive moving of tags


#' Simple interactive editing of tag clouds
#' 
#' A minimalistic editor for object of the tagcloud class.
#' 
#' \code{tagcloud} provides a minimalistic editor for tag clouds produced by
#' \code{\link{tagcloud}} function. After \code{editor.tagcloud} is called, the
#' tag cloud is plotted. First click selects the tag to be moved. The second
#' click sends the tag such that its left lower corner is at the position
#' indicated by the mouse. Right-clicking terminates the program.
#' 
#' @param boxes An object of the tagcloud class, returned by the
#' \code{\link{tagcloud}} function.
#' @return An object of the \code{tagcloud class} with the modified positions
#' of the tags.
#' @author January Weiner <january.weiner@@gmail.com>
#' @seealso \code{\link{tagcloud}}
#' @keywords tags tag clouds editing
#' @examples
#' 
#' \dontrun{
#' data( gambia )
#' terms <- gambia$Term
#' tagcloud( terms )
#' boxes <- editor.tagcloud( boxes )
#' 
#' }
#' 
#' @export editor.tagcloud
editor.tagcloud <- function( boxes ) {

  plot.tagcloud( boxes, with.box= T )
  nstep <- 10

  while( 1 ) {
    xvec <- as.vector( sapply( 1:nrow( boxes ), function( x ) seq( boxes[x,"x"], boxes[x,"x"] + boxes[x,"w"], length.out= nstep ) )  )
    yvec <- as.vector( sapply( 1:nrow( boxes ), function( x ) seq( boxes[x,"y"], boxes[x,"y"] + boxes[x,"h"], length.out= nstep ) )  )
    catf( "Please click on the label you want to move\n" )
    catf( "(right-click to finish)\n" )
    i <- identify( xvec, yvec, n= 1, plot= F )
    if( length( i ) == 0 ) break
    i <- as.integer( i / nstep ) + 1
    if( length( i ) == 0 ) break
    catf( "Please click on the new position for:\n" )
    catf( "%s\n", boxes$tags[i] )
    xy <- locator( 1 )
    debugpr( xy )
    boxes[i,"x"] <- xy[1]
    boxes[i,"y"] <- xy[2]
    plot( boxes, with.box= T )
  }
  plot( boxes )
  return( invisible( boxes ) )
}


