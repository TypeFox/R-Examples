`%but%` <- function( fun, x ){
  UseMethod( "%but%", x ) 
}

`%but%.default` <- function(fun,x){
  match.fun(fun)
}

`%but%.list` <- function( fun, x ){ 
  fun <- match.fun( fun )
  allArgs <- formals(fun)
  
  okArgs <- x[ names(x) %in% names(allArgs) ]
  allArgs[ names( okArgs ) ] <- okArgs  
  formals(fun) <- allArgs
  fun
}

`%but%.character` <- function( fun, x ){ 
  fun <- match.fun( fun )
  if( nchar(x) == 0) return(fun)
  
  allArgs <- formals(fun)
  test_is_logical <-  sapply( allArgs, is.logical )
  shortArgs <- substring( names(allArgs[test_is_logical]), 1,1)
  
  chars <- gregexpr( "[!-]?[a-zA-Z]", x )[[1]]
  chars <- substring( x, chars, chars + attr(chars, "match.length") - 1)
  for( current in chars ){
    actualChar <- gsub("[^a-zA-Z]", "", current )
    if( ! actualChar %in% shortArgs ) {
      warning( sprintf("No option with first letter `%s` in function", actualChar ) )
      next
    }
    
    ### if find a "!" ,  set the option to the opposite of the default
    allArgs[test_is_logical][[ which(actualChar == shortArgs)  ]] <- if( length(grep("!", current) ) ){
       !allArgs[test_is_logical][[ which(actualChar == shortArgs) ]]
    } else length(grep("-", current)) == 0
  }
  formals(fun) <- allArgs
  fun
} 


