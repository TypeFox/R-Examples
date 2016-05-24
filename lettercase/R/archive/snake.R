#' snake_case
#' 
#' Convert object to a snake_case object
#' 
#' @param x object to converted to snake case  
#' 
#' @examples
#'   snake( c("One Flew Over The Cuckoo's Nest", "Demon Box") )
#' 
#' @export
#' @rdname snake
#' 
#' @aliases snake_case

snake_case <- setClass( 'snake_case', contains='character' )

setMethod( 'initialize', 'snake_case', where=.GlobalEnv
  , function( .Object, x ) .Object@.Data <- .snake_case(x) 
)

snake_case( "one flew")


.snake_case <- function(x) { 

  if( ! is.character(x) ) stop( as.character(sys.call())[-1], ' is not character' )
  x <- tolower(x)
  x <- gsub( '\\s+', '_', x, perl=TRUE )
  
  class(x) <- "snake"
  x
  
}


#' @rdname snake
#' @export
#' @aliases as.snake
as.snake <- function(x) snake(x)


setOldClass( 'snake' )

setAs( 'character', 'snake', snake )


as( letters, 'snake' )
coerce(  )
#' character_to_snake <- 