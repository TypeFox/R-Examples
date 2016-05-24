# plots an object of the tagcloud class (or another suitable data frame)
#' @rdname tagcloud
#' @export
plot.tagcloud <- function( x, family= NULL, add= FALSE, with.box= FALSE, col= NULL, sel=NULL, ... ) {


  if( ! any( class( x ) == "data.frame" ) || ! any( class( x ) == "tagcloud" ) ) {
    stop( "x must be an object of class tagcloud or data.frame" )
  }

  boxes <- x

  # important: asp=1 guarantees that w and h of the boxes are
  # interchangeable
  if ( ! add ) {
    plot.new()
    old.par <- par( mar= c( 0, 0, 0, 0 ) )
    plot.window( xlim= c( 0, 1 ), ylim= c( 0, 1 ), asp= 1, ... )
  }


  if ( !missing( family ) ) {
    if ( length( family ) != nrow( x ) || length( family ) != 1 ) {
      stop( "Incorrect length of the family vector" )
    }
    boxes$family <- family
  }

  if ( !missing( col )) {
    if ( length( col ) != nrow( x ) || length( col ) != 1 ) {
      stop( "Incorrect length of the family vector" )
    }
    boxes$color  <- col
  }

  if ( !missing( sel ) ) {
    boxes <- boxes[sel,,drop=F]
  }
  
  for (i in 1:nrow(boxes)) {
    if (with.box)
      rect( boxes[i, "x" ], boxes[i, "y" ], 
      boxes[i, "x" ] + boxes[i, "w" ], boxes[i, "y" ] + boxes[i, "h" ] ) ;

    if ( boxes[i,"srt"] == 0 ) srt <- boxes[i,"vertical"] * 90
    else                       srt <- boxes[i,"srt" ]

    text( boxes[i,"x"] + boxes[i,"w"]/2, 
          boxes[i,"y"] + boxes[i,"h"]/2, 
          boxes$tags[i], cex= boxes[i,"cex"], family= boxes$family[i], 
          srt= srt, col= boxes$color[i] )
  }

  if( ! add ) par( old.par )
}


