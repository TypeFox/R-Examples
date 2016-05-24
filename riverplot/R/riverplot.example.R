#' Generate an example for riverplot
#'
#' The plotting functions in the riverplot package work on an object of the
#' riverplot class. This function returns an object of the riverplot class to
#' demonstrate how such an object (which is actually a simple list) can be
#' created.
#' 
#' 
#'
#' @title Generate an example for riverplot
#'
#'
#' @author January Weiner <january.weiner@@gmail.com>
#'
#' @export
#'
#' @examples
#' x <- riverplot.example()
#' plot( x )

riverplot.example <- function( ) {
  
  ret <- list( 

    nodes= data.frame(
      ID=LETTERS[1:8],
      x = c( 1, 2, 2, 3, 3, 4, 5, 1 ),
      labels= c( NA, NA, "Node C", rep( NA, 4 ), "Node H" ),
      stringsAsFactors= FALSE ),

  # nodes=  c(   B=2,   A=1,   Q=1, C=2,   D=3,   E=3, F=4, G=5 ),
   styles= list( 
     A=list( 
       col=  "#00990099",
       lty=0,
       textcol= "white"
       ),
     H=list(
       col= "#FF000099",
       textcol= "white"
     ),
     B= list( col= "#00006699", textcol= "white" ),
     F= list( col= "yellow" ),
     D= list( col= "#00FF0099" )
   )
  )

  ret$edges <- data.frame(
    N1=   c( "A", "A", "A", "H", "H", "H", "B", "B", "C", "C", "C" ),
    N2=   c( "B", "C", "D", "D", "F", "G", "D", "F", "D", "E", "F" ),
    Value= c( 10,  20,   5,  10,  10,  20,   5,  10,  20,  15,  10 ),
    stringsAsFactors= F )

  rownames( ret$nodes ) <- ret$nodes$ID
  class( ret ) <- c( class( ret), "riverplot" )
  return( ret )
}


