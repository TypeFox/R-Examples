#' http request helpers
#'
#' Various helpers as functions that allow chaining together options to
#' pass to e.g., \code{\link{Get}}
#'
#' @param .data A request object
#' @param seconds Number of seconds to wait for a response until giving up. Can
#' not be less than 1 ms (aka: 0.001 seconds).
#' @param agent An agent character string
#' @param user A user name
#' @param password A password
#' @param type Type of authenticaion. Default: basic
#' @param ... Further args
#'
#' @name helpers
NULL

#' @export
#' @rdname helpers
Progress <- function(.data){
  .data <- as.req(.data)
  .data <- modifyList(.data, list(config = c(progress())))
  .data
}

#' @export
#' @rdname helpers
Verbose <- function(.data){
  .data <- as.req(.data)
  .data <- modifyList(.data, list(config = c(verbose())))
  .data
}

#' @export
#' @rdname helpers
Timeout <- function(.data, seconds){
  .data <- as.req(.data)
  .data <- modifyList(.data, list(config = c(timeout(seconds = seconds) )))
  .data
}

#' @export
#' @rdname helpers
User_agent <- function(.data, agent){
  .data <- as.req(.data)
  .data <- modifyList(.data, list(config = c(user_agent(agent = agent))))
  .data
}

#' @export
#' @rdname helpers
Authenticate <- function(.data, user, password, type = "basic"){
  .data <- as.req(.data)
  .data <- modifyList(.data, list(config = c( authenticate(user = user, password = password, type = type) )))
  .data
}

#' @export
#' @rdname helpers
Query <- function(.data, ...){
  .data <- as.req(.data)
  args <- list(...)
  .data <- modifyList(.data, list(query = args))
  .data$parse <- TRUE
  .data
  # Get(.data)
}

#' @export
#' @rdname helpers
Body <- function(.data, ...){
  .data <- as.req(.data)
  args <- list(...)
  .data <- modifyList(.data, list(body = args))
  .data$parse <- TRUE
  Put(.data)
}

# query  <- function(.data=list(), ...){
#   query_(.data, .dots = lazyeval::lazy_dots(...))
# }
#
# query_  <- function(.data=list(), ..., .dots){
#   dots <- lazyeval::all_dots(.dots, ...)
#   args <- sapply(bb, "[[", "expr")
#   structure(args, class = "query")
# }
