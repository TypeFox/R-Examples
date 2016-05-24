predict.probit <- function ( object, newdata = NULL, type = "link", ... ) {
   
   if( is.null( newdata ) ) {
      mm <- model.matrix( object )
   } else {  # modified copy from predict.lm()
      tt <- terms( object )
      Terms <- delete.response( tt )
      m <- model.frame( Terms, newdata, xlev = object$xlevels )
      mm <- model.matrix( Terms, m )
   }
   result <- drop( mm %*% coef( object ) )
   if( type == "response" ) {
      result <- pnorm( result )
   } else if( type != "link" ) {
      stop( "argument 'type' must be either 'link' or 'response'" )
   }
   
   return( result )
}