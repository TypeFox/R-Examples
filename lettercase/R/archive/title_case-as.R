#' @include TitleCase.R
#' @name as
#' @rdname TitleCase
setAs( 'character', 'TitleCase', function(from) str_title_case(from) )

#' @name as
#' @rdname TitleCase

setAs( 'SnakeCase', 'TitleCase', function(from) str_title_case(from) )


#' @examples 
#'   as.title_case( "one flew over the cuckoo's_nest" )
#' @rdname TitleCase

as.title_case <- function(x) as(x, 'TitleCase' )       
