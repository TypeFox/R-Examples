which.article <- function(profiler)
{
  if (profiler == 'fpmpi' || profiler == 'mpip')
    return("An")
  else if (profiler == 'tau')
    return("A")
  else
    stop("Unknown profiler")
}



#' Printing
#' 
#' Print and show methods for \code{prof} class objects.
#' 
#' @docType methods
#' @param x,object
#' A \code{prof} class object
#' @param ...
#' extra arguments
#' 
#' @seealso \code{ \link{prof-class}, \link{read.prof} }
#' @keywords Methods
#' @name prof-print
#' @rdname prof-print
NULL



#' @rdname prof-print
#' @export
setMethod("print", signature(x="prof"),
  function(x, ...)
  {
    cat(sprintf(paste(which.article(x@profiler), x@profiler, "profiler object:\n")))
    print(x@parsed)
  }
)



#' @rdname prof-print
#' @export
setMethod("show", signature(object="prof"),
  function(object){
    print(object)
  }
)

