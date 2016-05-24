#' Assign with type checking
#'
#' @export
#' @name typed-assign
#' @param x symbol
#' @param ... additional conditions taking the following forms:
#'
#'     1. \code{fun = v}, i.e. \code{fun(x)} must be equal \code{v}.
#'
#'     2. \code{cond}, i.e. \code{cond(x)} must be \code{TRUE}.
#'
#'     3. a \code{function} like \code{function(x) mean(x) <= 5.0}
#' @param value value to be assigned
#' @examples
#' \dontrun{
#' x <- 10L
#' atomic(x) <- 20
#' numeric(x) <- 10
#' numeric(x, length = 10L) <- 1:10
#'
#' cond1 <- function(x) mean(x) <= 5
#' numeric(x, cond1) <- 0:9
#' }
`atomic<-` <- function(x, ..., value) {
  check(x, value, is.atomic, ...)
}

#' @rdname typed-assign
#' @export
`integer<-` <- function(x, ..., value) {
  check(x, value, is.integer, ...)
}

#' @rdname typed-assign
#' @export
`numeric<-` <- function(x, ..., value) {
  check(x, value, is.numeric, ...)
}

#' @rdname typed-assign
#' @export
`double<-` <- function(x, ..., value) {
  check(x, value, is.double, ...)
}

#' @rdname typed-assign
#' @export
`logical<-` <- function(x, ..., value) {
  check(x, value, is.logical, ...)
}

#' @rdname typed-assign
#' @export
`character<-` <- function(x, ..., value) {
  check(x, value, is.character, ...)
}

#' @rdname typed-assign
#' @export
`raw<-` <- function(x, ..., value) {
  check(x, value, is.raw, ...)
}

#' @rdname typed-assign
#' @export
`complex<-` <- function(x, ..., value) {
  check(x, value, is.complex, ...)
}

#' @rdname typed-assign
#' @export
`matrix<-` <- function(x, ..., value) {
  check(x, value, is.matrix, ...)
}

#' @rdname typed-assign
#' @export
`array<-` <- function(x, ..., value) {
  check(x, value, is.array, ...)
}

#' @rdname typed-assign
#' @export
`list<-` <- function(x, ..., value) {
  check(x, value, is.list, ...)
}

#' @rdname typed-assign
#' @export
`pairlist<-` <- function(x, ..., value) {
  check(x, value, is.pairlist, ...)
}

#' @rdname typed-assign
#' @export
`envir<-` <- function(x, ..., value) {
  check(x, value, is.environment, ...)
}

#' @rdname typed-assign
#' @export
`name<-` <- function(x, ..., value) {
  check(x, value, is.name, ...)
}

#' @rdname typed-assign
#' @export
`symbol<-` <- function(x, ..., value) {
  check(x, value, is.symbol, ...)
}

#' @rdname typed-assign
#' @export
`call<-` <- function(x, ..., value) {
  check(x, value, is.call, ...)
}

#' @rdname typed-assign
#' @export
`factor<-` <- function(x, ..., value) {
  check(x, value, is.factor, ...)
}

#' @rdname typed-assign
#' @export
`fun<-` <- function(x, ..., value) {
  check(x, value, is.function, ...)
}

#' @rdname typed-assign
#' @export
`expression<-` <- function(x, ..., value) {
  check(x, value, is.expression, ...)
}

#' @rdname typed-assign
#' @export
`language<-` <- function(x, ..., value) {
  check(x, value, is.language, ...)
}

#' @rdname typed-assign
#' @export
`object<-` <- function(x, ..., value) {
  check(x, value, is.object, ...)
}

#' @rdname typed-assign
#' @export
`table<-` <- function(x, ..., value) {
  check(x, value, is.table, ...)
}

#' @rdname typed-assign
#' @export
`recursive<-` <- function(x, ..., value) {
  check(x, value, is.recursive, ...)
}

#' @rdname typed-assign
#' @export
`vector<-` <- function(x, ..., value) {
  check(x, value, is.vector, ...)
}

#' @rdname typed-assign
#' @export
`data.frame<-` <- function(x, ..., value) {
  check(x, value, is.data.frame, ...)
}

#' @rdname typed-assign
#' @export
`null<-` <- function(x, ..., value) {
  check(x, value, is.null, ...)
}

#' @rdname typed-assign
#' @export
`check<-` <- function(x, ..., value) {
  check(x, value, type = , ...)
}
