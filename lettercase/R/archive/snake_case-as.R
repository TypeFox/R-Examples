#' @include SnakeCase.R
#' @rdname SnakeCase
#' @name as

setAs( 'character', 'SnakeCase', function(from) str_snake_case(from) )

#' @rdname SnakeCase
#' @name as

setAs( 'TitleCase', 'SnakeCase', function(from) str_snake_case(from) )


# @examples 
#   as.snake_case( "one flew over the cuckoo's_nest" )
#' @rdname SnakeCase
#' @export

as.snake_case <- function(x) as(x, 'SnakeCase' )       
