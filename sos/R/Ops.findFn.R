
#' binary operators for RSiteSearch objects
Ops.findFn <- function( e1, e2 ){
  if( nargs() != 2){
    stop( sprintf( "unary operator '%s' not implemented", .Generic) )
  }

# making sure both e1 and e2 are RSiteSearch objects
  rework <- function( x ){
    if( !inherits( x, "findFn" ) ){
      if( inherits( x, "character" ) && length(x) ){
        findFn( paste( x, collapse= "+" ), quiet = TRUE )
      } else{
        msg <- paste("findFn objects can only be combined",
                     "with either findFn objects or",
                     "character vectors" )
        stop(msg)
      }
    } else{
      x
    }
  }
  e1 <- rework( e1 )
  e2 <- rework( e2 )
  e12 <- switch( .Generic,
         "&" = intersectFindFn( e1, e2 ),
         "|" = unionFindFn( e1, e2),
         stop( sprintf(
           "operator '%s' not implemented for findFn objects",
                       .Generic)  )
         )
  e12
}

