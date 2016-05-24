"mksub" <-
function(x, start, end, id) {
  orig <- as.character(substitute(x))[[1]]
  a <- list()  # make a copy of the column attributes
  for (n in names(x)) {
    tmp <- attributes(x[[n]])  # strip out unneeded attributes
    srp <- names(tmp) %in% c("dim", "dimnames", "names")
    a[[n]] <- if (any(srp))
      tmp[!srp]
    else
      tmp
  }
  if (missing(id))
    id <- NULL
  if ("id" %in% names(x)) {
    x$id <- as.character(x$id)
    if (is.null(id)) {
      id <- unique(x$id)
      if (length(id) > 1) {
        id <- id[1]
        warning(gettextf("more than one unique id found in %s; using %s",
                         sQuote(sprintf("%s$id", orig)), sQuote(id)))
      }
    }
    if (!is.null(id)) {
      x <- x[x$id == id,]
      x$id <- NULL
      attr(x, "id") <- as.character(id)[[1]]
      attr(x, "name") <- c(getstnname(id), orig)[[1]]
    }
    if (nrow(x) <= 0)
      stop(gettextf("no data at %s", sQuote(sprintf("id=%s", id))))
  } else {
    if (!is.null(id))
      warning(gettextf("%s must have an %s column to select an ID",
                       sQuote(orig), sQuote("id")))
    attr(x, "name") <- orig
  }
  sc <- seas.df.check(x, orig)  # check to see if it is a seasonal df
  if (missing(start))
    start <- NULL
  if (inherits(start, c("POSIXct", "Date"))) {
    date.range <- range(x$date)
    if (!is.null(start) && start < date.range[1])
      warning(gettextf("%s is before the minimum date in %s of %s",
                       sQuote("start"), sQuote(orig), sQuote(date.range[1])))
    if (missing(end)) {
      if (is.null(start))
        end <- as.Date("0000-01-01")
      else
        end <- as.Date(format(start, "%Y-12-31"))
    }
    if (end > date.range[2])
      warning(gettextf("%s is after the maximum date in %s of %s",
                       sQuote("end"), sQuote(orig), sQuote(date.range[2])))
    date <- x$date
  } else {  # start and end are years
    date <- as.integer(format(x$date, "%Y"))
    date.range <- range(date)
    if (!is.null(start) && start < date.range[1])
      warning(gettextf("%s is before the minimum year in %s of %s",
                       sQuote("start"), sQuote(orig), sQuote(date.range[1])))
    if (missing(end))
      end <- 0
    else if (end > date.range[2])
      warning(gettextf("%s is after the maximum year in %s of %s",
                       sQuote("end"), sQuote(orig), sQuote(date.range[2])))
  }
  if (is.null(start)) {
    start <- min(date)
    end <- max(date)
  }
  if (start > end)  # one year
    x <- x[date == start,]
  else  # range of years
    x <- x[start <= date & date <= end,]
  for (n in names(x))  # copy attributes back
    attributes(x[[n]]) <- a[[n]]
  x
}
