#' Declare symbols
#' @name declare
#' @export
#' @param ... Symbols to declare
#' @param .envir \code{environment} to store the symbols
#' @examples
#' declare(x,y=numeric(),z=integer())
declare <- function(...,.envir=parent.frame()) {
  args <- dots(...)
  Map(function(name,value) {
    if(nchar(name) == 0L) {
      name <- deparse(value)
      value <- NULL
    }
    assign(name, eval(value, .envir), envir = .envir)
  },getnames(args),args)
  invisible(NULL)
}
