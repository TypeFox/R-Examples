##' Compute frequencies
##'
##' @param x factor
##' @param useNA useNA
##' @param propNA propNA
##' @param cum logical
##' @param addmargins addmargins
##' @author David Hajage
##' @keywords internal
freq <- function(x, useNA = c("no", "ifany", "always"), propNA = TRUE, cum = FALSE, addmargins = FALSE) {
  rnames <- as.character(as.list(substitute(list(x)))[-1])
  n <- n.table(x, useNA = useNA, margin = 0, addmargins = addmargins)
  p <- p.table(x, useNA = useNA, propNA = propNA, margin = 0, addmargins = FALSE)[[1]]
  if (addmargins)
    p <- as.table(c(p, Total = 1))
  if (cum) {
    n.cum <- cumsum(n)
    p.cum <- cumsum(p)
  } else {
    n.cum <- NULL
    p.cum <- NULL
  }
  results <- cbind(n, p, n.cum, p.cum)
  attr(results, "lgroup") <- list(names(n), rnames)
  attr(results, "n.lgroup") <- list(1, nrow(results))
  class(results) <- c("freq", "matrix")
  results
}

##' Compute frequencies (data.frame input)
##'
##' @importFrom Hmisc label
##' @param df data.frame
##' @param useNA useNA
##' @param propNA propNA
##' @param cum logical
##' @author David Hajage
##' @keywords internal
freq.data.frame <- function(df, useNA = c("no", "ifany", "always"), propNA = TRUE, cum = FALSE, addmargins = FALSE, label = FALSE) {
  if (cum)
    addmargins <- FALSE
  dfl <- as.list(df)

  if (!label)
    rnames <- names(dfl)
  else
    rnames <- sapply(dfl, Hmisc:::label.default)
  
  results <- lapply(dfl, freq, useNA = useNA, propNA = propNA, cum = cum, addmargins = addmargins)
  nrows <- sapply(results, nrow)
  results <- rbind.list(results)

  attr(results, "n.lgroup") <- list(1, nrows)
  attr(results, "lgroup") <- list(unlist(lapply(df, function(x) {
    if (useNA[1] == "always" | ((useNA[1] == "ifany") & any(is.na(x))))
      lev <- levels(addNA(x))
    else 
      lev <- levels(x)
    if (addmargins)
      lev <- c(lev, "Total")
    lev})), rnames)
  class(results) <- c("freq", "matrix")
  attr(results, "df") <- df
  results
}

##' Ascii for freq object.
##'
##' Ascii method for freq object (internal).
##'
##' @export
##' @method ascii freq
##' @import ascii
##' @param x a freq object
##' @param format see \code{?ascii} in \code{ascii} package
##' @param digits see \code{?ascii} in \code{ascii} package
##' @param include.rownames see \code{?ascii} in \code{ascii} package
##' @param rownames see \code{?ascii} in \code{ascii} package
##' @param include.colnames see \code{?ascii} in \code{ascii} package
##' @param header see \code{?ascii} in \code{ascii} package
##' @param lgroup see \code{?ascii} in \code{ascii} package
##' @param n.lgroup see \code{?ascii} in \code{ascii} package
##' @param ... other arguments passed to \code{ascii}
##' @author David Hajage
##' @keywords univar
ascii.freq <- function(x, format = "nice", digits = 3, include.rownames = FALSE, include.colnames = TRUE, header = TRUE, lgroup = attr(x, "lgroup"), n.lgroup = attr(x, "n.lgroup"), ...) {
  class(x) <- class(x)[-1]
  ascii(x, include.colnames = include.colnames, include.rownames = include.rownames, header = header, lgroup = lgroup, n.lgroup = n.lgroup, format = format, digits = digits, ...)
}

##' Print freq object.
##'
##' Print freq object (internal).
##'
##' @export
##' @import ascii
##' @method print freq
##' @param x a freq object
##' @param type type of output (see \code{?ascii} in \code{ascii}
##' package)
##' @param lstyle see \code{?ascii} in \code{ascii} package
##' @param ... other arguments passed to \code{ascii}
##' @author David Hajage
##' @keywords univar
print.freq <- function(x, type = "rest", lstyle = "", ...) {
  print(ascii.freq(x, lstyle = lstyle, ...), type = type)
  invisible(x)
}

##' as.data.frame for freq object.
##'
##' as.data.frame for freq object (internal).
##'
##' @export
##' @param x a freq object
##' @param ... not used
##' @author David Hajage
##' @keywords internal
as.data.frame.freq <- function(x, ...) {
  xx <- unclass(x)
  var <- unlist(mapply(rep, attr(x, "lgroup")[[2]], attr(x, "n.lgroup")[[2]], SIMPLIFY = FALSE))
  levels <- attr(x, "lgroup")[[1]]
  
  data.frame(var = var, levels = levels, xx, row.names = NULL, check.names = FALSE)
}

##' Test if \code{x} is an freq object
##'
##' @export
##' @param x a freq object
is.freq <- function(x)
  inherits(x, "freq")
