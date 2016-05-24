
# prints a summary for the tagcloud class
#' @rdname tagcloud
#' @export
summary.tagcloud <- function( object, ... ) {

  boxes <- object

  ratio <- boxes.ratio( boxes )
  bb    <- c( min( boxes[,"x"] ), min( boxes[,"y"] ),
              max( boxes[,"x"] + boxes[,"w"] ), max( boxes[,"y"] + boxes[,"h"] ) )
   
  ret <-  list( 
      n= nrow( boxes ), coverage= ratio,
      weights.range= range( boxes$weights ),
      bb= bb, algorithm= attr( boxes, "algorithm" ),
      scale= attr( boxes, "scale" )
    ) 

  class( ret ) <- c( "tagcloudsummary", class( ret ) )
  return( ret ) 
}

print.tagcloud <- function( x, ... ) 
  print.tagcloudsummary( summary.tagcloud( x ), ... ) 

print.tagcloudsummary <- function( x, ... ) {

  catf( "Object of the tagcloud class\n" )
  catf( "  number of tags: %d\n", x$n )
  catf( "  value range: %.2f - %.2f\n", x$weights.range[1], x$weights.range[2] )
  catf( "  coverage: %.0f %%\n", 100 * x$coverage )
  catf( "  bb: (%.2f, %.2f), (%.2f, %.2f)\n", x$bb[1], x$bb[2], x$bb[3], x$bb[4] )
  catf( "  width: %.2f height: %.2f\n", x$bb[3] - x$bb[1], x$bb[4] - x$bb[2] )
  catf( "  algorithm: %s\n", x$algorithm )
  catf( "  scale: %.3f\n", x$scale )
 
}
