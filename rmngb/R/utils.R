asInteger <- function(x)
    as.integer(round(x))

invDiag <- function(x) {
    vec <- integer(x ^ 2)
    vec[cumsum(rep(x - 1, x)) + 1] <- 1
    matrix(vec, nrow = x)
}

"%out%" <- function(x, y)
    ! x %in% y

nameToString <- function(x)
    deparse(substitute(x))

rmAttr <- function(x, except = "class") {
    nAttr <- names(attributes(x))
    attributes(x) <- attributes(x)[except]
    x
}

Rnumber <- function()
    paste((v <- R.Version())$major, v$minor, sep = ".")

Rname <- function()
    R.Version()$nickname

closest <- function(x, y, tol = +Inf)
    (d <- abs(x - y)) == min(d) & (d < tol)

locf <- function(x)
    (x2 <- na.omit(x))[length(x2)]

colClasses <- function(x)
    unlist(lapply(lapply(x, class), function(x) x[length(x)]))

nDistinct <- function(x)
    length(table(x))

interleave <- function(x, y) {
  iX <- 2 * seq_along(x) - 1
  iY <- 2 * seq_along(y)
  c(x, y)[order(c(iX, iY))]
}
