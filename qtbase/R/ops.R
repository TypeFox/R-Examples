### Override R operators for Qt objects

Ops.RQtObject <- function(e1, e2) {
  fun <- paste("operator", .Generic, sep = "")
  ans <- try({
    if (missing(e2))
      qinvokeStatic(Qt$QGlobalSpace, fun, e1)
    else qinvokeStatic(Qt$QGlobalSpace, fun, e1, e2)
  }, silent = TRUE)
  if (inherits(ans, "try-error"))
    ans <- NextMethod()
  ans
}

"[.RQtObject" <- function(x, i, j, ... , drop = TRUE) {
  if (length(list(...)) || !missing(j) || !missing(drop))
    stop("all arguments except 'x' and 'i' are ignored")
  qinvoke("operator[]", x, i)
}

"%<<%" <- function(x, y) UseMethod("%<<%")
"%<<%.RQtObject" <- function(x, y) qinvoke("operator<<", x, y)

"%>>%" <- function(x, y) UseMethod("%>>%")
"%>>%.RQtObject" <- function(x, y) qinvoke("operator>>", x, y)
