##' Compute summary statistics according to a factor
##'
##' @param x numeric object
##' @param by factor
##' @param funs functions
##' @param ... passed to funs
##' @param useNA useNA
##' @param addmargins addmargins
##' @author David Hajage
##' @keywords internal
summarize.by <- function(x, by, funs = c(mean, sd, quantile, n, na), ..., useNA = c("no", "ifany", "always"), addmargins = FALSE) {
  if (!is.numeric(x))
    stop("x doit etre numerique")  

  if (useNA[1] != "no") {
    by <- addNA(by)
  }
  
  if (!is.character(funs)) {
    funs <- as.character(as.list(substitute(funs)))
    funs <- funs[funs != "c" & funs != "list"]
  }

  fun <- function(df, funs = funs, ...) summarize(df[, 1], funs = funs, ...)
  df <- data.frame(x, by, check.names = FALSE)
  names(df) <- c(names(df)[1], "by")
  results <- list(dlply(df, .(by), fun, funs = funs, ...))

  lgroup <- lapply(results, function(x) {
    ## y <- names(x)
    y <- as.character(attr(x, "split_labels")[, 1])
    y[is.na(y)] <- "NA"
    y
  })
  n.lgroup <- NULL

  if (addmargins) {
    results[[1]] <- c(results[[1]], list(summarize(x, funs = funs, ...)))
    lgroup <- c(lgroup[[1]], "Total")
  }
  
  class(results) <- c("summarize.by", "list")
  attr(results, "split_type") <- NULL
#  attr(results, "split_labels") <- NULL
  attr(results, "lgroup") <- lgroup
  attr(results, "n.lgroup") <- n.lgroup
  attr(results, "revert") <- FALSE
  
  results
}

##' Compute summary statistics according to a factor (data.frame input)
##'
##' @importFrom Hmisc label
##' @import ascii
##' @param df data.frame
##' @param by data.frame
##' @param funs fuctions
##' @param ... passed to funs
##' @param useNA useNA
##' @param addmargins addmargins
##' @param revert whether to regroup factors or numeric variables when crossing factor with numeric variables
##' @author David Hajage
##' @keywords internal
summarize.data.frame.by <- function(df, by, funs = c(mean, sd, quantile, n, na), ..., useNA = c("no", "ifany", "always"), addmargins = FALSE, revert = FALSE, test = FALSE, test.summarize = test.summarize.auto, show.test = display.test, plim = 4, show.method = TRUE, label = FALSE) {
  if (!is.character(funs)) {
    funs <- as.character(as.list(substitute(funs)))
    funs <- funs[funs != "c" & funs != "list"]
  }

  if (useNA[1] != "no") {
    for (i in 1:ncol(by))
      by[, i] <- addNA(by[, i])
  }

  if (!label) {
    namescoldf <- names(df)
    namescolby <- names(by)
  } else {
    namescoldf <- sapply(df, Hmisc:::label.default)
    namescolby <- sapply(by, Hmisc:::label.default)
  }
  
  nby <- ncol(by)
  namesby <- paste("by", 1:nby, sep = ".")
  fun <- function(df, funs = funs, ...) summarize.data.frame(df[, 1:(ncol(df)-nby), drop = FALSE], funs = funs, ...)
  dff <- data.frame(df, by, check.names = FALSE)
  names(dff) <- c(names(dff)[1:(length(names(dff))-nby)], namesby)
  results <- NULL
  for (i in 1:nby) {
    dfi <- dff
    if (useNA[1] == "no") {
      dfi <- dfi[!is.na(by[, i]), ]
    }
    results <- c(results, list(dlply(dfi, namesby[i], fun, funs = funs, ...)))
  }
  names(results) <- colnames(by)

  r <- lapply(results, function(x) {
    ## y <- names(x)
    y <- as.character(attr(x, "split_labels")[, 1])    
    y[is.na(y)] <- "NA"
    if (addmargins)
      y <- c(y, "Total")
    y
  })

  if (addmargins) {
    tot.results <- summarize.data.frame(df, funs = funs, ...)
    results <- lapply(results, function(x) c(x, list(tot.results)))
  }
  
  if (revert) {
    lgroup <- list(rep(namescoldf, length(unlist(r))), unlist(r), namescolby)
    nr <- nrow(results[[1]][[1]])
    n.lgroup <- list(1, nr, nr*sapply(r, length))
    xx <- NULL
    for (i in 1:length(results)) {
      for (j in 1:length(results[[i]])) {
        xx <- rbind(xx, results[[i]][[j]])
      }
    }
    results <- xx

    rgroup <- NULL
    n.rgroup <- NULL
    if (test) {
      tests <- NULL
      for (i in 1:nby) {
        dfi <- df
        byi <- by
        if (useNA[1] == "no") {
          dfi <- dfi[!is.na(by[, i]), , drop = FALSE]
          byi <- byi[!is.na(by[, i]), , drop = FALSE]
        }
        tests <- c(tests, list(sapply(split(dfi, byi[, i]), function(x) {
          d <- stack(x)
          show.test(test.summarize(d$values, d$ind), digits = plim, method = show.method)
        })))
      }
      names(tests) <- colnames(by)
      rgroup <- unlist(tests)
      n.rgroup <- rep(ncol(df), length(rgroup))
    }
  } else {
    lgroup <- list(unlist(lapply(r, rep, ncol(df))), rep(namescoldf, ncol(by)), namescolby)
    nr <- rep(sapply(results, length), each = ncol(df))
    nrr <- sapply(results, length)*ncol(df)
    n.lgroup <- list(1, nr, nrr)
    results <- do.call(rbind, lapply(results, function(x) do.call(interleave.matrix, x)))

    rgroup <- NULL
    n.rgroup <- NULL
    if (test) {
      tests <- NULL
      for (i in 1:nby) {
        dfi <- df
        byi <- by
        if (useNA[1] == "no") {
          dfi <- dfi[!is.na(by[, i]), , drop = FALSE]
          byi <- byi[!is.na(by[, i]), , drop = FALSE]        
        }
        tests <- c(tests, list(sapply(dfi, function(x) {
          show.test(test.summarize(x, byi[, i]), digits = plim, method = show.method)
        })))
      }
      names(tests) <- colnames(by)
      rgroup <- unlist(tests)
      n.rgroup <- n.lgroup[[2]]
    }
  }

  class(results) <- c("summarize.by", "list")
  attr(results, "split_type") <- NULL
  attr(results, "split_labels") <- NULL
  attr(results, "lgroup") <- lgroup
  attr(results, "n.lgroup") <- n.lgroup
  attr(results, "rgroup") <- rgroup
  attr(results, "n.rgroup") <- n.rgroup
  attr(results, "revert") <- revert

  attr(results, "df") <- df
  attr(results, "by") <- by
  
  results
}

##' Ascii for summarize.by object.
##'
##' Ascii method for summarize.by object (internal).
##'
##' @export
##' @method ascii summarize.by
##' @import ascii
##' @param x a summarize.by object
##' @param format see \code{?ascii} in \code{ascii} package
##' @param digits see \code{?ascii} in \code{ascii} package
##' @param include.rownames see \code{?ascii} in \code{ascii} package
##' @param include.colnames see \code{?ascii} in \code{ascii} package
##' @param header see \code{?ascii} in \code{ascii} package
##' @param lgroup see \code{?ascii} in \code{ascii} package
##' @param n.lgroup see \code{?ascii} in \code{ascii} package
##' @param rgroup see \code{?ascii} in \code{ascii} package
##' @param n.rgroup see \code{?ascii} in \code{ascii} package
##' @param rstyle see \code{?ascii} in \code{ascii} package
##' @param ... other arguments passed to \code{ascii}
##' @author David Hajage
##' @keywords univar
ascii.summarize.by <- function(x, format = "nice", digits = 5, include.rownames = FALSE, include.colnames = TRUE, header = TRUE, lgroup = attr(x, "lgroup"), n.lgroup = attr(x, "n.lgroup"), rgroup = attr(x, "rgroup"), n.rgroup = attr(x, "n.rgroup"), rstyle = "d", ...) {
  ascii(unclass(x), lgroup = lgroup, n.lgroup = n.lgroup, rgroup = rgroup, n.rgroup = n.rgroup, rstyle = rstyle, include.rownames = include.rownames, include.colnames = include.colnames, header = header, format = format, digits = digits, ...)
}

##' Print summarize.by object.
##'
##' Print summarize.by object (internal).
##'
##' @export
##' @method print summarize.by
##' @import ascii
##' @param x a summarize.by object
##' @param type type of output (see \code{?ascii} in \code{ascii} package)
##' @param ... other arguments passed to \code{ascii}
##' @author David Hajage
##' @keywords internal
print.summarize.by <- function(x, type = "rest", lstyle = "", ...) {
  print(ascii.summarize.by(x, lstyle = lstyle, ...), type = type)
  ## invisible(x)
}

##' as.data.frame for summarize.by object.
##'
##' as.data.frame for summarize.by object (internal).
##'
##' @export
##' @param x a summarize.by object
##' @param ... not used
##' @author David Hajage
##' @keywords internal
as.data.frame.summarize.by <- function(x, ...) {
  xx <- unclass(x)
  if (!attr(x, "revert")) {
    lgroup <- attr(x, "lgroup")
    n.lgroup <- attr(x, "n.lgroup")
    lgroup[[2]] <- unlist(mapply(rep, lgroup[[2]], each = n.lgroup[[2]], SIMPLIFY = FALSE))
    lgroup[[3]] <- unlist(mapply(rep, lgroup[[3]], n.lgroup[[3]], SIMPLIFY = FALSE))
    xx <- data.frame(by = lgroup[[3]], var = lgroup[[2]], levels = lgroup[[1]], xx, row.names = NULL, check.names = FALSE)
  } else {
    lgroup <- attr(x, "lgroup")
    n.lgroup <- attr(x, "n.lgroup")
    lgroup[[2]] <- unlist(mapply(rep, lgroup[[2]], each = n.lgroup[[2]], SIMPLIFY = FALSE))
    lgroup[[3]] <- unlist(mapply(rep, lgroup[[3]], n.lgroup[[3]], SIMPLIFY = FALSE))
    xx <- data.frame(by = lgroup[[3]], levels = lgroup[[2]], var = lgroup[[1]], xx, row.names = NULL, check.names = FALSE)
  }
  xx
}

##' Test if \code{x} is a summarize.by object
##'
##' @param x a summarize.by object
##' @author David Hajage
##' @keywords internal
is.summarize.by <- function(x)
  inherits(x, "summarize.by")
