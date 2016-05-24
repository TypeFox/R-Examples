#' @title Is class
#' @description Tests if an object is of a given class
#' @param x an R object
#' @param name string giving the class to be tested
#' @export
#' @examples
#' is_class("test_me", "character") # TRUE
#' is_class(1:10, "numeric") # TRUE
#' 
#' y = 'hello'
#' class(y) = "hello"
#' is_class(y, 'hello')
is_class <- function(x, name=NULL) {
  if (is.null(name))
    stop("\n'name' is missing with no default")
  if (is_string(name)) {
    if (class(x) == name) TRUE else FALSE    
  } else {
    stop("\n'name' must be a string")
  }
}
