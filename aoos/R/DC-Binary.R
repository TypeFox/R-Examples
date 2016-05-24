#' Binary-class
#' 
#' This is a virtual class to be contained in other class definitions. It can be used to define binary operators, e.g. \code{+} or \code{-}, inside an aoos class definition (\code{\link{defineClass}}).
#' 
#' @exportClass Binary
#' 
#' @details At the moment you can define binary operators as methods by naming them as \code{.<binaryOperator>} (see the example). This is implemented for the following operators: \code{+, -, *, /, \%\%, ^, <, >, ==, >=, <=, &}.
#' @rdname Binary
#' 
#' @examples
#' Rational <- defineClass("Rational", contains = c("Show", "Binary"), {
#'
#'   numer <- 0
#'   denom <- 1
#'   .g <- 1
#'   
#'   .gcd <- function(a, b) if(b == 0) a else Recall(b, a %% b)
#'
#'   init <- function(numer, denom) {
#'     .self$.g <- .gcd(numer, denom)
#'     .self$numer <- numer / .g
#'     .self$denom <- denom / .g
#'   }
#'   
#'   show <- function() {
#'     cat(paste0(.self$numer, "/", .self$denom, "\n"))
#'   }
#'   
#'   ".+" <- function(that) {
#'     Rational(numer = numer * that$denom + that$numer * denom,
#'              denom = denom * that$denom)
#'   }
#'   
#'   neg <- function() {
#'     Rational(numer = -.self$numer,
#'              denom = .self$denom)
#'   }
#'   
#'   ".-" <- function(that) {
#'      .self + that$neg()
#'   }
#'   
#' })
#'
#' rational <- Rational(2, 3)
#' rational + rational
#' rational$neg()
#' rational - rational
setClass("Binary", contains = "VIRTUAL")

binaryMethodNames <- function() {
  list(c("+", ".+"),
       c("-", ".-"),
       c("*", ".*"),
       c("/", "./"),
       c("%%", ".%%"),
       c("^", ".^"),
       c("<", ".<"),
       c(">", ".>"),
       c("==", ".=="),
       c(">=", ".>="),
       c("<=", ".<="),
       c("&", ".&"))
}

makeBinaryMethod <- function(funName) {
  force(funName)
  function(e1, e2) {
    e <- as.environment(e1)
    if (exists(funName, envir = e)) {
      if (missing(e2)) get(funName, envir = e)() 
      else get(funName, envir = e)(e2)
    } else {
      NextMethod()
    }
  }
}

for (methodNamePair in binaryMethodNames()) {
  setMethod(methodNamePair[1], c("Binary", "ANY"), makeBinaryMethod(methodNamePair[2]))
}
