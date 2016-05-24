#' Suppresses all output except for errors.
#' 
#' Evaluates an expression and suppresses all output except for errors, 
#' meaning: prints, messages, warnings and package startup messages.
#' 
#' @param expr [valid R expression]\cr
#'   Expression. 
#' @return Return value of expression invisibly. 
#' @export 
#' @examples
#' suppressAll({
#'   print("foo")
#'   message("foo")
#'   warning("foo")
#' })
suppressAll = function(expr) {
  capture.output({
    z = suppressWarnings(
      suppressMessages(
        suppressPackageStartupMessages(force(expr))                  
      )
    )
  })
  invisible(z)
}
