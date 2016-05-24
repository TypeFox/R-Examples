#' Upstreamness
#' 
#' @name upstream
#' @param x an object of class "decompr" as created using the load_tables_vectors() function from the decompr package.
#' @export
#' @import decompr
#' @examples 
#' # load the decompr package
#' library(decompr)
#' 
#' # load example data
#' data(leather)
#' 
#' # create a leontief decomposed data set
#' l <- load_tables_vectors(inter,
#'                          final,
#'                          countries,
#'                          industries,
#'                          out        )
#'  
#'  # apply upstream
#'  upstream( l )

upstream <- function ( x ) {
  
  solve( diag(x$GN) - (t(t(x$A)*x$X)/x$X) ) %*% matrix(1, nrow=x$GN)
  
}
