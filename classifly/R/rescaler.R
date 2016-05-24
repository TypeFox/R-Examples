rescaler <- function(x, type="sd", ...) UseMethod("rescaler", x)

#' @export
rescaler.default <- function(x, type="sd", ...) {
  switch(type,
    rank = rank(x, ...),
    var = ,
    sd = (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE),
    robust = (x - median(x, na.rm=TRUE)) / mad(x, na.rm=TRUE),
    I = x,
    range = (x - min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE))
  )
}

#' @export
rescaler.data.frame <- function(x, type="sd", ...) {
  continuous <- sapply(x, is.numeric)
  x[continuous] <- lapply(x[continuous], rescaler, type=type, ...)
  x
}

#' @export
rescaler.matrix <- function(x, type="sd", ...) {
  apply(x, 2, rescaler, type=type, ...)
}
