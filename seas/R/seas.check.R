"seas.df.check" <-
function(x, orig, var=NULL) {
  if (missing(orig))
    orig <- as.character(substitute(x))[[1]]
  if (!inherits(x, "data.frame"))
    stop(gettextf("%s is not a %s object", sQuote(orig), sQuote("data.frame")))
  if (!"date" %in% names(x))
    stop(gettextf("a %s column must exist in %s", sQuote("date"), sQuote(orig)))
  if (!inherits(x$date, c("POSIXct", "Date")))
    stop(gettextf("%s$date must be either %s or %s class",
                  sQuote(orig), sQuote("Date"), sQuote("POSIXct")))
  sc <- list()
  sc$id <- attr(x, "id")
  if (is.null(sc$id) && !is.null(x$id)) {
    sc$id <- unique(as.character(x$id))
    if (length(sc$id) > 1)
      stop(gettextf("more than one unique id found in %s",
                    sQuote(sprintf("%s$id", orig))))
  }
  if ("name" %in% names(attributes(x)))
    sc$name <- attr(x, "name")
  if (is.null(sc$name) &&!is.null(sc$id))
    sc$name <- getstnname(sc$id)
  if (is.null(sc$name))
    sc$name <- orig
  sc$year.range <- as.integer(format(range(x$date, na.rm=TRUE), "%Y"))
  sc$calendar <- if ("date" %in% names(x) &&
                    "calendar" %in% names(attributes(x$date)))
    attr(x$date, "calendar")
  else if ("calendar" %in% names(attributes(x)))
    attr(x, "calendar")
  else NULL
  sc$start.day <- if ("date" %in% names(x) &&
                      "start.day" %in% names(attributes(x$date)))
    attr(x$date, "start.day")
  else if ("start.day" %in% names(attributes(x)))
    attr(x, "start.day")
  else NULL
  sc$main <- .seastitle(id=sc$id, name=sc$name, orig=orig, range=sc$year.range)
  if (!is.null(var)) {
    vn <- var %in% names(x)
    if (!all(vn))
      stop(
        sprintf(
          ngettext(sum(!vn), "%s is not found in %s", "%s are not found in %s"),
          paste(sprintf("%s", sQuote(var[!vn])), collapse=", "), sQuote(orig)))
    vn <- sapply(x[,var, drop=FALSE],
                 function(x)sum(is.finite(x))) < 1
    if (any(vn))
      stop(
        sprintf(
          ngettext(sum(vn), "%s from %s has no data", "%s from %s have no data"),
          paste(sprintf("%s", sQuote(var[vn])), collapse=", "), sQuote(orig)))
    var <- var[1]
    sc$units <- attr(x[[var]], "units")
    sc$long.name <- attr(x[[var]], "long.name")
    sc$ylab <- .seasylab(var, sc$long.name, sc$units)
  }
  invisible(sc)
}

"seas.sum.check" <-
function(x, orig, var, norm, year.filter, ann.only) {
  if (is.null(getOption("seas.main")))
  if (missing(orig))
    orig <- as.character(substitute(x))[[1]]
  if (!inherits(x, "seas.sum"))
    stop(gettextf("%s is not a %s object", sQuote(orig), sQuote("seas.sum")))
  if (missing(var)) {
    var <-x$prime[[1]]
  }
  x$var <- var
  if (!(var %in% x$var))
    stop(gettextf("%s not found in %s", sQuote(var),sQuote(orig)))
  if (missing(ann.only))
    ann.only <- FALSE
  if (!ann.only) {
    if (inherits(norm, "matrix")) {
      if (!all(dim(norm) == dim(x$seas[1:2])))
        stop(gettextf("%s does not have the same dimensions as %s",
                      sQuote(norm), sQuote(sprintf("%s$days", orig))))
      x$days <- norm
    } else if (inherits(norm, "array")) {
      if (!all(dim(norm) == dim(x$seas)))
        stop(gettextf("%s does not have the same dimensions as %s",
                      sQuote(norm), sQuote(sprintf("%s$seas", orig))))
      x$norm <- norm
    } else {
      if (norm == "active") {
        if (!"active" %in% names(x))
          stop(gettextf("%s not found in %s", sQuote("active"),sQuote(orig)))
        x$norm <- x$active
      } else if (!"norm" %in% names(x)) {
        x$norm <- array(x$days, dim=dim(x$seas), dimnames=dimnames(x$seas))
      }
    }
  }
  if (!missing(year.filter)) {
    s <- x$years %in% year.filter
    x$ann <- x$ann[s,, drop=FALSE]
    if (!ann.only) {
      x$seas <- x$seas[s,,, drop=FALSE]
      x$norm <- x$norm[s,,, drop=FALSE]
      if (x$a.cut)
        x$active <- x$active[s,,, drop=FALSE]
    }
    x$na <- x$na[s,, drop=FALSE]
    x$days <- x$days[s,, drop=FALSE]
    x$years <- x$years[s]
  }
  if (length(x$years) < 1)
    stop("no data in %s", sQuote(orig))
  invisible(x)
}
