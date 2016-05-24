##This file contains some small helper functions.

makePval<-function(zs) {
  2*pnorm(abs(zs),lower.tail=FALSE)
}

##' Returns \code{formula} attribute of an \code{xbal} object.
##'
##' @param x An \code{xbal} object.
##' @param ... Ignored.
##' @return The formula corresponding to \code{xbal}.
##' @export
formula.xbal<-function(x,...){
  attr(x,"fmla")
}

# ============================================================================
# = withOptions helper provides a safe way to temporarily override options() =
# ============================================================================
withOptions <- function(optionsToChange, fun) {
  oldOpts <- options()
  do.call(options, optionsToChange)
  tryCatch(fun(), finally = do.call(options, oldOpts))
}

##Our own version of these to handle the signif stars.
###print.ftable<-function (x, digits = getOption("digits"), ...) {
###  write.ftable(x, quote = FALSE, digits = digits)
###}
###
###write.ftable<-function (x, file = "", quote = TRUE, append = FALSE, digits = getOption("digits"),justify.labels="right",justify.data="right",...)
###{
###    r <- RItools:::format.ftable(x, quote = quote, digits = digits,justify.labels=justify.labels,justify.data=justify.data,...)
###    cat(t(r), file = file, append = append, sep = c(rep(" ",
###        ncol(r) - 1), "\n"))
###    invisible(x)
###}
###
###format.ftable<-function (x, quote = TRUE, digits = getOption("digits"), justify.labels="left",justify.data="right", ...)
###{
###    if (!inherits(x, "ftable"))
###        stop("'x' must be an \"ftable\" object")
###    charQuote <- function(s) if (quote)
###        paste("\"", s, "\"", sep = "")
###    else s
###    makeLabels <- function(lst) {
###        lens <- sapply(lst, length)
###        cplensU <- c(1, cumprod(lens))
###        cplensD <- rev(c(1, cumprod(rev(lens))))
###        y <- NULL
###        for (i in rev(seq_along(lst))) {
###            ind <- 1 + seq.int(from = 0, to = lens[i] - 1) *
###                cplensD[i + 1]
###            tmp <- character(length = cplensD[i])
###            tmp[ind] <- charQuote(lst[[i]])
###            y <- cbind(rep(tmp, times = cplensU[i]), y)
###        }
###        y
###    }
###    makeNames <- function(x) {
###        nmx <- names(x)
###        if (is.null(nmx))
###            nmx <- rep("", length.out = length(x))
###        nmx
###    }
###    xrv <- attr(x, "row.vars")
###    xcv <- attr(x, "col.vars")
###    LABS <- cbind(rbind(matrix("", nrow = length(xcv), ncol = length(xrv)),
###        charQuote(makeNames(xrv)), makeLabels(xrv)), c(charQuote(makeNames(xcv)),
###        rep("", times = nrow(x) + 1)))
###    DATA <- rbind(if (length(xcv))
###        t(makeLabels(xcv)), rep("", times = ncol(x)), format(unclass(x),
###        digits = digits))
###    cbind(apply(LABS, 2, format, justify = justify.labels), apply(DATA,
###        2, format, justify = justify.data))
###}

##' Select variables, strata, and statistics from a \code{xbal} object
##'
##' If any of the arguments are not specified, all the of relevant
##' items are included.
##'
##' @param x The \code{xbal} object, the result of a call to
##'   \code{\link{xBalance}}
##' @param vars The variable names to select.
##' @param strata The strata names to select.
##' @param stats The names of the variable level statistics to select.
##' @param tests The names of the group level tests to select.
##' @param ... Other arguments (ignored)
##'
##' @return A \code{xbal} object with just the appropriate items
##'   selected.
##'
##' @export
subset.xbal <- function(x,
                        vars   = NULL,
                        strata = NULL,
                        stats  = NULL,
                        tests  = NULL,
                        ...) {

  res.dmns <- dimnames(x$results)
  ov.dmns <- dimnames(x$overall)

  if (is.null(strata)) {
    strata <- res.dmns$strata
  }

  if (is.null(vars)) {
    vars <- res.dmns$vars
  }

  if (is.null(stats)) {
    stats <- res.dmns$stat
  }

  if (is.null(tests)) {
    tests <- ov.dmns$tests
  }

  res <- x$results[vars, stats, strata, drop = F]
  ovr <- x$overall[strata, tests, drop = F]

  attr(res, "originals") <- attr(x$results, "originals")[res.dmns$vars %in% vars]

  tmp <- list(results = res, overall = ovr)
  class(tmp) <- c("xbal", "list")

  return(tmp)
}
