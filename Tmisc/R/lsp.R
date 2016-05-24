#' List objects in package
#' 
#' Lists functions and how to call them for any package.
#' 
#' @author Stephen Turner
#' @keywords keywords
#' 
#' @param package The name of the package you're examining.
#' @param ... further arguments to be passed to \code{lsf.str}.
#' 
#' @return A list of functions and how to call them for any package.
#' @import utils
#' 
#' @examples
#' \dontrun{
#' lsp(Tmisc, pattern="un")
#' }
#' 
#' @export
lsp <-function(package, ...) {
    package <- deparse(substitute(package))
    lsf.str(paste0("package:", package), ...)
}
