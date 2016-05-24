checkFuncArgs <- function( func, checkArgs, argName, funcName ) {
   ## is the 'func' a function?
   if( !is.function( func ) ) {
      stop( "argument '", argName, "' of function '", funcName,
         "' is not a function" )
   }

   funcArgs <- names( formals( func ) )

   if( length( funcArgs ) > 1 ) {
      a <- charmatch( funcArgs[ -1 ], checkArgs )
      if( sum( !is.na( a ) ) == 1 ) {
         stop( "argument '", funcArgs[ -1 ][ !is.na( a ) ],
            "' of the function specified in argument '", argName,
            "' of function '", funcName,
            "' (partially) matches the argument names of function '",
            funcName, "'. Please change the name of this argument" )
      } else if( sum( !is.na( a ) ) > 1 ) {
         stop( "arguments '",
            paste( funcArgs[ -1 ][ !is.na( a ) ], collapse = "', '" ),
            "' of the function specified in argument '", argName,
            "' of function '", funcName,
            "' (partially) match the argument names of function '",
            funcName, "'. Please change the names of these arguments" )
      }
   }
   return( NULL )
}

