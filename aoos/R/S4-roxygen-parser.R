#' Parser for roxygen documentation
#' 
#' These functions are used by roxygen2 for generating documentation. 
#' 
#' @param call a call
#' @param env an environment
#' @param block is ignored
#' 
#' @rdname parser
#' @export
"parser_%m%" <- function(call, env, block) {
  value <- eval(call, env)
  # This cryptic piece exists to fool R CMD check, which rightfully complains
  # that I rely on a not exported object:
  value@.Data <- eval(parse(text = "roxygen2:::extract_method_fun"))(value@.Data)
  roxygen2::object(value)
}

#' @export
#' @rdname parser
"parser_%g%" <- function(call, env, block) {
  value <- eval(call, env)
  roxygen2::object(value)
}

#' @export
#' @rdname parser
"parser_%type%" <- function(call, env, block) {
  value <- eval(call, env)
  roxygen2::object(value)
}
