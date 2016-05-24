#' Create S3 methods
#' 
#' Helper function to create S3 method that preserves (most) 
# attributes (including class).
#' 
#' @param generic,object arguments passed to \code{\link{NextMethod}}
#' 
#' @keywords internal
#' 
#' @return 
#' Function with arguments "x", "i", and "...". 
#' The "i" argument can be thought of as "index" if the method defined is 
#' for example the "["-function.  
create_s3_method <- function(generic = NULL, object = NULL){
  function(x, i, ...) {
    r <- NextMethod(generic = generic, object = object)
    mostattributes(r) <- attributes(x)
    r
  }
}