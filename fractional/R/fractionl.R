#' Representation of a numeric vector in vulgar fractional form
#'
#' The object is flagged so that if it is coerced to \code{character},
#' or printed, the numerical quantities are represented by a rational
#' approximation.  In other respects the numerical object behaves as
#' normally.
#'
#' @param x A numeric object
#' @param eps An absolute error tolerance
#' @param maxConv An upper limit on the number of convergents to use in the
#' continued fractions.
#' @param sync A logical value.  Should the numerical value be changed to match
#' the rational approximation, as closely as possible with floating point, (TRUE)? Or,
#' should it be left and used in its original state (FALSE)?
#' @param ... Currently ignored.
#'
#' @return A numeric object of class \code{"fractional"}.
#' @seealso \code{\link[MASS]{fractions}} for a similar functionality.
#' @export
#'
#' @examples
#' (M <- solve(cbind(1, contr.helmert(5))))
#' (Mf <- fractional(M))     ## print method right justifies
#' (Mc <- as.character(Mf))  ## print method left justifies
#' (Mn <- numerical(Mc))
#' set.seed(123)
#' u <- matrix(runif(10), 2, 5)
#' (uf <- fractional(u))
#' (us <- fractional(u, sync = TRUE))  ## may look different!
#' unfractional(uf) - unfractional(us)  ## rational approximation errors
fractional <- function(x, eps = 1.0e-6, maxConv = 20,
                       sync = FALSE) {
  if(sync) {
    ax <- attributes(x)
    nd <- rat(as.vector(x), eps = eps, maxConv = maxConv)
    x <- nd[, "Pn"]/nd[, "Qn"]
    eps <- .Machine$double.eps
    attributes(x) <- ax
  }
  structure(x, eps = eps, maxConv = maxConv, sync = sync,
            class = c("fractional", class(x)))
}

#' Demote a fractional object back to a numeric one
#'
#' Given an object of class \code{"fractional"} this simple function
#' removes the attributes that signal that it is to be treated as a
#' fractional object, thus returning it to its original \code{numeric}
#' status alone
#'
#' @param x A \code{"fractional"} object
#'
#' @return A simple \code{numeric} object like \code{x}
#' @export
#'
#' @examples
#' (tst <- fractional(matrix(0:9/10, 2, 5)))
#' (tst <- unfractional(tst))
unfractional <- function(x) {
  x <- unclass(x)
  attr(x, "eps") <- attr(x, "maxConv") <- attr(x, "sync") <- NULL
  x
}

#' @describeIn ratr Workhorse function for a single value
#' @export
.ratr <- function(x, eps = 1.0e-6, maxConv = 20) {
  PQ1 <- c(1, 0)
  PQ2 <- c(floor(x), 1)
  r <- x - PQ2[1]
  i <- 0
  while((i <- i+1) < maxConv && abs(x - PQ2[1]/PQ2[2]) > eps) {
    b <- floor(1/r)
    r <- 1/r - b
    PQ0 <- PQ1
    PQ1 <- PQ2
    PQ2 <- b*PQ1 + PQ0
  }
  return(c(PQ2, i-1))
}

#' Calculate Rational Approximation Using Continued Fraction Methods
#'
#' This is a behind-the-scenes function not likely to be used other than
#' internally within the package.  It computes the rational approximations for
#' each value in the principal argument.
#'
#' @param x A numeric vector for which rational approximations are required.
#'
#' @param eps An absolute error tolerance on the approximation
#'
#' @param maxConv An upper limit on the number of convergents that the continued
#'  fraction expansion may employ.  The fraction is terminated once the desired
#'  accuracy is met (or the upper limit is about to be exceeded).
#'
#' @return A 3 column matrix giving, respectively, the numerators, denominators
#' and number of convergents needed to achieve the error tolerance, in the columns
#'
#' @examples
#' fractional(base::pi)
#' ratr(base::pi)
#'
#' set.seed(123)
#' (u <- matrix(runif(10), 2, 5))
#' (ru <- ratr(u, eps = 1.0e-3, maxConv = 6))
#' (abs_error <- matrix(abs(u - ru[, 1]/ru[, 2]), 2, 5))
#'
#' @seealso \code{\link{rat}} which has the same functionality, but is coded in \code{C++}.
#' @export
ratr <- function(x, eps = 1.0e-6, maxConv = 20) {
  structure(t(sapply(x, .ratr, eps = eps, maxConv = maxConv)),
            dimnames = list(NULL, c("Pn", "Qn", "n")))
}

getAttr <- function(x) UseMethod("getAttr")
getAttr.default <- function(x) list(eps = 1.0e-6, maxConv = 20, sync = FALSE)
getAttr.fractional <- function(x) attributes(x)[c("eps", "maxConv", "sync")]


#' Method for the group generic function for the arithmetic operators
#'
#' Provides arithmetic operations for numeric objects or of class \code{"fractional"}.
#'
#'
#' @param e1 A numeric objet, possibly of class \code{"fractional"}
#' @param e2 A numeric objet, possibly of class \code{"fractional"}
#'
#' @return The result of the arithmetic operation, flagged as class
#' \code{"fractional"}
#' @export
#'
#' @examples
#' (M <- fractional(1:10/7))
#' M + 1
#' 1 + M + M^2
#'
Ops.fractional <- function (e1, e2) {
  ax <- getAttr(e1)
  e1 <- unclass(e1)
  if (!missing(e2)) {
    ax2 <- getAttr(e2)
    ax <- list(eps = max(.Machine$double.eps, min(ax$eps, ax2$eps)),
               maxConv = max(ax$maxConv, ax2$maxConv),
               sync = ax$sync & ax2$sync)
    e2 <- unclass(e2)
  }
  res <- NextMethod(.Generic)
  if(typeof(res) == "logical") {
    res
  } else {
    with(ax, fractional(res, eps = eps, maxConv = maxConv,
                        sync = sync))
  }
}

#' Method for the group generic function for the elementary mathematical functions
#'
#' Allows graceful operations with mathematical functions.
#'
#' @param x  A numerical object flagged as \code{fractional}
#' @param ... Passed on to further methods (but usually not required)
#'
#' @return A numeric object with the results of the computations, but NOT flagged
#' as of class \code{"fractional"}.
#' @export
#'
#' @examples
#' (M <- fractional(solve(cbind(1, contr.helmert(5)))))
#' (M0 <- abs(M)*sign(M))  ## fractional flag lost
Math.fractional <- function(x, ...) {
  x <- unclass(x)
  attr(x, "eps") <- attr(x, "maxConv") <- attr(x, "sync") <- NULL
  NextMethod(.Generic, x, ...)
}

#' @describeIn fractional S3 method for coercion to character,
#' producing an object inheriting from class \code{"charFrac"}
#'
#' @export
as.character.fractional <- function (x, eps = attr(x, "eps"),
                                     maxConv = attr(x, "maxConv"), ...) {
  x <- unclass(x)
  ax <- attributes(x)
  rx <- rat(as.vector(x), eps = eps, maxConv = maxConv)
  fractions <- sub("/1$", "", paste(rx[, "Pn"], rx[, "Qn"], sep = "/"))
  ax$maxConv <- ax$eps <- ax$sync <- NULL
  attributes(fractions) <- ax
  class(fractions) <- "charFrac"
  fractions
}

#' @describeIn fractional Print method for class \code{"charFrac"} objects, unquoted.
#' @export
print.charFrac <- function(x, ...) {
  y <- x
  x <- gsub("^0$", ".", unclass(x))
  NextMethod("print", quote = FALSE, ...)
  invisible(y)
}

#' Convert a fractional object to the equivalent numeric object
#'
#' Convert an object of class \code{"fractional"} or \code{"charFrac"} to a purely
#' numeric object.  This is effectively a method function for the \code{.Primitive}
#' generic function \code{as.numeric} but written as a separate function for purely
#' technical reasons.
#'
#' @param vulgar character string form of a class 'fractional' object.
#' @return A \code{numeric} object as represented by its (usually \code{fractional}) display.
#' @export
#' @examples
#' suppressPackageStartupMessages(library(dplyr))
#' m <- 2*diag(5)
#' m[abs(row(m) - col(m)) == 1] <- -1
#' m  ## How much roundoff error does inverting entail?
#' (mi  <- solve(m) %>% fractional) ## patterned inverse
#' mi * max(denominators(mi)) ## clearer pattern
#' m1 <- solve(mi)
#' range(m1 - m)           ## roundoff still present
#' m2 <- m1 %>% numerical  ## remove roundoff error - hopefully!
#' identical(m2, m)        ## no roundoff
numerical <- function(vulgar) {
  UseMethod("numerical")
}

#' @describeIn numerical Method for \code{"fractional"} objects
#' @export
numerical.fractional <- function(vulgar) {
  numerators(vulgar)/denominators(vulgar)
}

#' @describeIn numerical Method for \code{"charFrac"} objects
#' @export
numerical.charFrac <- function(vulgar) {
  numerators(vulgar)/denominators(vulgar)
}

#' @describeIn numerical Default method for \code{numerical} generic
#' @export
numerical.default <- function(vulgar) {
  nd <- rat(vulgar)
  vulgar[] <- nd[, "Pn"]/nd[, "Qn"]
  vulgar
}

#' @describeIn fractional Print method for \code{"fractional"} objects
#' @export
print.fractional <- function (x, ...) {
  x0 <- x
  y <- gsub("^0$", ".", as.character.fractional(x))
  y <- format(y, justify = "right")
  ax <- attributes(x)
  ax$class <- ax$eps <- ax$maxConv <- ax$sync <- NULL
  x <- do.call("structure", c(list(y), ax))
  NextMethod("print", quote = FALSE, ...)
  invisible(x0)
}

#' @export
`[.fractional` <- function(x, ...) {
  attr <- attributes(x)
  r <- unclass(x)[...]
  attributes(r) <-
    c(attributes(r), attr[!(names(attr) %in%
                              c("dim", "dimnames", "names"))])
  r
}

#' @export
`[<-.fractional` <- function(x, ..., value) {
  x <- unclass(x)
  r <- NextMethod()
  with(attributes(x), fractional(r, eps, maxConv, sync))
}

## Housekeeping when fractional is detached and unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("fractional", libpath)
}
