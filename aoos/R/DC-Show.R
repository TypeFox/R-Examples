#' Show class
#' 
#' This is a virtual class to be contained in other class definitions. It overrides the default show method and is intended to be used with the aoos class system (\code{\link{defineClass}}). The show method will simply look for a method \code{show} defined as member of a class definition.
#' 
#' 
#' @seealso \link{defineClass}
#' @export
#' @rdname Show
#' @examples
#' ClassWithShowMethod <- defineClass("ClassWithShowMethod", contains = "Show", {
#'   show <- function() print(summary(.self))
#' })
#' 
#' ClassWithShowMethod()
setClass("Show", contains = "VIRTUAL")

#' @export
#' @rdname Show
#' 
#' @param object an object inheriting from \code{Show}
setMethod("show", "Show", function(object) {
  object$show()
})
