#' @name rmall
#' @keywords remove all
#' @author Sven E. Templer
#' @title Remove All Objects from Global Environment
#' @description 
#' Remove all objects from the global environment.
#' @param ... Arguments forwarded to \code{ls} to get all objects.
#' @seealso
#' \link{rm}, \link{ls}
#' @examples
#' #
#' 
#' a <- b <- letters
#' ls()
#' rmall()
#' ls()
#' 
#' #

#' @export rmall
rmall <- function (...) {
  n <- ls(envir=.GlobalEnv, ...)
  rm(list=n, envir=.GlobalEnv)
  invisible(NULL)
}
