"position<-" <- function(x, value) {
  x <- as.ordered(x)
  if (length(levels(x)) != length(value))
    stop("length(levels(x)) != length(value)")
  if (!is.positioned(x))
    if.R(r=class(x) <- c("positioned", class(x)),
         s=oldClass(x) <- c("positioned", oldClass(x)))
  attr(x, "position") <- value
  x
}


position <- function(x) {
  if(is.positioned(x))
    attr(x, "position")
  else
    if (inherits(x,  "factor"))
      as.numeric(seq(along=levels(x)))
    else
      as.numeric(x)
}

unpositioned <- function(x, ...) {
  if (!is.positioned(x)) return(x)
  oldClass(x) <- oldClass(x)["positioned" != oldClass(x)]
  attr(x, "position") <- NULL
  x
}

is.numeric.positioned <- function(x, ...) { ## S-Plus
  ## S-Plus requires this and it doesn't get dispatched
  ## because is.numeric is not a method.
  ## It must be called explicitly
  if (is.positioned(x)) FALSE
  else
    is.numeric(x)
}

as.position <- function(x) {
  if (is.numeric(x))
    x
  else
    position(x)[x]
}

"[.positioned" <- function (x, ..., drop = FALSE) {
  if.R(s={drop=FALSE}, r={}) ## S-Plus has an argument matching problem
  position.x <- position(x)
  x <- unpositioned(x)
  y <- NextMethod("[")
  if (drop) position(y) <- position.x[!is.na(match(levels(x), levels(y)))]
  else  position(y) <- position.x
  y
}

as.positioned <- function(x) {
  if (is.positioned(x)) x
  else {
    z <- as.ordered(x)
    position(z) <- position(z)
    z
  }
}

is.positioned <- function(x)
  inherits(x, "positioned")

positioned <- function(x, ..., value) {
  x <- ordered(x, ...)
  position(x) <-
    if (!missing(value)) value
    else position(x)
  x
}

print.positioned <- function(x, ...) {
  xx <- x
  attr(x, "position") <- NULL
  oldClass(x) <- oldClass(x)[-1]
  NextMethod("print")
  cat("position:", paste(position(xx), collapse = " < "), "\n")
  invisible(xx)
}

unique.positioned <- function(x, incomparables = FALSE, ...) {
  y <- unique(unpositioned(x), incomparables = FALSE, ...)
  position(y) <- position(x)
  y
}

if.R(s=
     setOldClass(c("positioned", "ordered", "factor"))
     ,
     r={})

as.numeric.positioned <- function(x, ...){
  if (is.numeric.positioned(x))
    x
  else
    position(x)[as.numeric(unpositioned(x))]
}

is.na.positioned <- function(x) {
  ## S-Plus requires this
  is.na(unpositioned(x))
}

## source("~/HH-R.package/HH/R/position.R")
