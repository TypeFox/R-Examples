#' get raw data
#'
#' @param .obj Data
#' @param .data A request object
#' @param .dots xxx
#' @param parse (logical) To parse or not
#' @param ... further args
#'
#' @name raw
NULL

#' @export
#' @rdname raw
raw <- function(.obj = list(), ...){
  raw_(.obj, .dots = lazyeval::lazy_dots(...), parse = TRUE)
}

#' @export
#' @rdname raw
raw_ <- function(.data = list(), ..., .dots, parse = TRUE){
  dots <- lazyeval::all_dots(.dots, ...)
  Get(.data, dots)
}
