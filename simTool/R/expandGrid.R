#' Creates a \code{data.frame} from All Combinations
#'
#'  Actually a wrapper for \code{\link{expand.grid}}, but
#'  character vectors will stay as characters.
#'
#'
#'@param \dots  vectors, factors or a list containing these.
#'@return See \code{\link{expand.grid}}
#'@author  Marsel Scheer
#'@seealso  \code{\link{expand.grid}}
#'@examples
#'
#'expandGrid(fun="rnorm", mean=1:4, sd=2:5)
#'
#'@export
expandGrid <-
function(...){
  expand.grid(..., stringsAsFactors=FALSE)
}
