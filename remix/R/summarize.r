##' Compute summary statistics
##'
##' @param x numeric object
##' @param funs functions
##' @param ... passed to funs
##' @author David Hajage
##' @keywords internal
summarize <- function(x, funs = c(mean, sd, quantile, n, na), ...) {
  if (!is.numeric(x))
    stop("x doit etre numerique")  

  if (!is.character(funs)) {
    funs <- as.character(as.list(substitute(funs)))
    funs <- funs[funs != "c" & funs != "list"]
  }

  if (length(funs) > 1)
    fun <- do.call(funs2fun, as.list(funs))
  else
    fun <- match.fun(funs)

  results <- fun(x, ...)
  class(results) <- c("summarize", "matrix")
  attr(results, "x") <- x
  results
}

##' Compute summary statistics (data.frame input)
##'
##' @importFrom Hmisc label
##' @param df a data.frame
##' @param funs functions
##' @param ... passed to funs
##' @author David Hajage
##' @keywords internal
summarize.data.frame <- function(df, funs = c(mean, sd, quantile, n, na), label = FALSE, ...) {
  if (!is.character(funs)) {
    funs <- as.character(as.list(substitute(funs)))
    funs <- funs[funs != "c" & funs != "list"]
  }

  dfl <- as.list(df)
  results <- t(sapply(dfl, summarize, funs = funs, ...))  
  class(results) <- c("summarize", "matrix")
  if (length(funs) == 1) {
    if (length(match.fun(funs)(1:10, ...)) == 1) {
      dim(results) <- rev(dim(results))
      ## rownames(results) <- names(dfl)
      colnames(results) <- funs
    }
  }
  if (!label)
    rownames(results) <- names(dfl)
  else
    rownames(results) <- sapply(dfl, Hmisc:::label.default)
  attr(results, "df") <- df
  results
}

##' Ascii for summarize object.
##'
##' Ascii method for summarize object (internal).
##'
##' @export
##' @method ascii summarize
##' @import ascii
##' @param x a summarize object
##' @param format see \code{?ascii} in \code{ascii} package
##' @param digits see \code{?ascii} in \code{ascii} package
##' @param include.rownames see \code{?ascii} in \code{ascii} package
##' @param include.colnames see \code{?ascii} in \code{ascii} package
##' @param header see \code{?ascii} in \code{ascii} package
##' @param ... other arguments passed to \code{ascii}
##' @author David Hajage
##' @keywords univar
ascii.summarize <- function(x, format = "nice", digits = 5, include.rownames = TRUE, include.colnames = TRUE, header = TRUE, ...) {
  if (is.null(nrow(x))) {
    x <- t(x)
  }
  class(x) <- class(x)[-1]
  ascii(x, include.rownames = include.rownames, include.colnames = include.colnames, header = header, format = format, digits = digits, ...)
}

##' Print summarize object.
##'
##' Print summarize object (internal).
##'
##' @export
##' @method print summarize
##' @import ascii 
##' @param x a summarize object
##' @param type type of output (see \code{?ascii} in \code{ascii}
##' package)
##' @param lstyle see \code{?ascii} in \code{ascii} package
##' @param ... other arguments passed to \code{ascii}
##' @author David Hajage
##' @keywords internal
print.summarize <- function(x, type = "rest", lstyle = "", ...) {
  print(ascii.summarize(x, lstyle = lstyle, ...), type = type)
  ## invisible(x)
}

##' as.data.frame for summarize object.
##'
##' as.data.frame for summarize object (internal).
##'
##' @export
##' @param x a summarize object
##' @param ... not used
##' @author David Hajage
##' @keywords internal
as.data.frame.summarize <- function(x, ...) {
  as.data.frame(unclass(x))
}

##' Test if \code{x} is a summarize object
##'
##' @param x a summarize object
##' @author David Hajage
##' @keywords internal
is.summarize <- function(x)
  inherits(x, "summarize")
