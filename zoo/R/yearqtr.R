## class creation
yearqtr <- function(x) structure(floor(4*x + .001)/4, class = "yearqtr")

## coercion to yearqtr: always go via numeric
as.yearqtr <- function(x, ...) UseMethod("as.yearqtr")
as.yearqtr.default <- function(x, ...) as.yearqtr(as.numeric(x))
as.yearqtr.numeric <- function(x, ...) structure(floor(4*x + .0001)/4, class = "yearqtr")
as.yearqtr.integer <- function(x, ...) structure(x, class = "yearqtr")

# as.jul.yearqtr <- function(x, ...) jul(as.Date(x, ...)) # jul is from tis
as.yearqtr.mondate <-
as.yearqtr.jul <- # jul is in tis package
as.yearqtr.timeDate <-
as.yearqtr.dates <-
as.yearqtr.Date <- 
as.yearqtr.POSIXt <- function(x, ...) as.yearqtr(as.yearmon(x))
as.yearqtr.yearqtr <- function(x, ...) x

as.yearqtr.factor <- function(x, ...) as.yearqtr(as.character(x), ...)
as.yearqtr.character <- function(x, format, ...) {
    non.na <- x[!is.na(x)]
    if (length(non.na) == 0) 
        return(structure(rep(NA, length(x)), class = "yearqtr"))
    if (missing(format) || format == "") {
        format <- if (all(regexpr("q", non.na) > 0))  { "%Y q%q"
        } else if (all(regexpr("Q", non.na) > 0)) { "%Y Q%q"
        } else "%Y-%q"
    }
    y <- if (regexpr("%[qQ]", format) > 0) {
        format <- sub("%q", "%m", format)
        y <- as.numeric(as.yearmon(x, format))
        m0 <- round(12 * (y %% 1))
        floor(y) + ifelse(m0 > 3, NA, m0/4)
    } else as.yearmon(x, format)
    as.yearqtr(y)
}
as.yearqtr.ti <- function(x, ...) as.yearqtr(as.Date(x), ...)

## coercion from yearqtr
# returned Date is the fraction of the way through the period given by frac
as.Date.yearqtr <- function(x, frac = 0, ...) {
  x <- unclass(x)
  if(all(is.na(x))) return(as.Date(x))
  year <- floor(x + .001)
  ix <- !is.na(year)
  month <- floor(12 * (x - year) + 1 + .5 + .001)
  dd.start <- as.Date(rep(NA, length(year)))
  dd.start[ix] <- as.Date(paste(year[ix], month[ix], 1, sep = "-")) 
  dd.end <- dd.start + 100 - as.numeric(format(dd.start + 100, "%d")) 
  as.Date((1-frac) * as.numeric(dd.start) + frac * as.numeric(dd.end), origin = "1970-01-01")
}
as.POSIXct.yearqtr <- function(x, tz = "", ...) as.POSIXct(as.Date(x), tz = tz, ...)
as.POSIXlt.yearqtr <- function(x, tz = "", ...) as.POSIXlt(as.Date(x), tz = tz, ...)
as.list.yearqtr <- function(x, ...) lapply(seq_along(x), function(i) x[i])
as.numeric.yearqtr <- function(x, ...) unclass(x)
as.character.yearqtr <- function(x, ...) format.yearqtr(x, ...)
as.data.frame.yearqtr <- function(x, row.names = NULL, optional = FALSE, ...) 
{
  nrows <- length(x)
  nm <- paste(deparse(substitute(x), width.cutoff = 500), collapse = " ")
  if (is.null(row.names)) {
    if (nrows == 0) 
        row.names <- character(0)
    else if(length(row.names <- names(x)) == nrows && !any(duplicated(row.names))) {
    }
    else if(optional) row.names <- character(nrows)
    else row.names <- seq_len(nrows)
  }
  names(x) <- NULL
  value <- list(x)
  if(!optional) names(value) <- nm
  attr(value, "row.names") <- row.names
  class(value) <- "data.frame"
  value
}


## other methods for class yearqtr
c.yearqtr <- function(...) {
    as.yearqtr(do.call("c", lapply(list(...), as.numeric)))
}

cycle.yearqtr <- function(x, ...) round(4 * (as.numeric(x) %% 1)) + 1

format.yearqtr <- function(x, format = "%Y Q%q", ...) 
{
    if (length(x) == 0) return(character(0))
	# like gsub but replacement and x may be vectors the same length
	gsub.vec <- function(pattern, replacement, x, ...) {
		y <- x
		for(i in seq_along(x)) {
			y[i] <- gsub(pattern, replacement[i], x[i], ...)
		}
		y
	}
	x <- as.yearqtr(x)
	x <- unclass(x)
	year <- floor(x + .001)
	qtr <- floor(4*(x - year) + 1 + .5 + .001)
    if (format == "%Y Q%q") return(paste(year, " Q", qtr, sep = ""))
    # TODO: speed up the following
	xx <- gsub.vec("%q", qtr, rep(format, length(qtr)))
	xx <- gsub.vec("%Y", year, xx)
	xx <- gsub.vec("%y", sprintf("%02d", as.integer(year %% 100)), xx)
	xx <- gsub.vec("%C", year %/% 100, xx)
	names(xx) <- names(x)
	xx
}


months.yearqtr <- function(x, abbreviate = FALSE) {
    months(as.Date(x), abbreviate = abbreviate)
}

quarters.yearqtr <- function(x, abbreviate = FALSE) {
    quarters(as.Date(x), abbreviate = abbreviate)
}


print.yearqtr <- function(x, ...) { 
    print(format(x), ...)
    invisible(x) 
}

"[.yearqtr" <- function (x, ..., drop = TRUE) 
{
    cl <- oldClass(x)
    class(x) <- NULL
    val <- NextMethod("[")
    class(val) <- cl
    val
}

"[[.yearqtr" <- function (x, ..., drop = TRUE) 
{
    cl <- oldClass(x)
    class(x) <- NULL
    val <- NextMethod("[[")
    class(val) <- cl
    val
}

MATCH.yearqtr <- function(x, table, nomatch = NA, ...)
    match(floor(4*as.numeric(x) + .001), floor(4*as.numeric(table) + .001), nomatch = nomatch, ...)

Ops.yearqtr <- function(e1, e2) {
    e1 <- as.numeric(as.yearqtr(e1))
    e2 <- as.numeric(as.yearqtr(e2))
    rval <- NextMethod(.Generic)
    if(is.numeric(rval)) rval <- yearqtr(rval)
    return(rval)
}


"-.yearqtr" <- function (e1, e2) 
{
    if (!inherits(e1, "yearqtr")) 
        stop("Can only subtract from yearqtr objects")
    if (nargs() == 1) 
	return(- as.numeric(e1))
    if (inherits(e2, "yearqtr")) 
        return(as.numeric(e1) - as.numeric(e2))
    if (!is.null(attr(e2, "class"))) 
      stop("can only subtract yearqtr objects and numbers from yearqtr objects")
    yearqtr(unclass(e1) - e2)
}

is.numeric.yearqtr <- function(x) FALSE

Axis.yearqtr <- function(x = NULL, at = NULL, ..., side, labels = NULL)
    axis.yearqtr(x = x, at = at, ..., side = side, labels = TRUE)


axis.yearqtr <- function (side, x, at, format, labels = TRUE, ..., N1 = 25, N2 = 7) {
    # If years in range > N1 then only years shown.  
    # If years in range > N2 then quarter ticks are not labelled.
    mat <- missing(at) || is.null(at)
    if (!mat) # at not missing
        x <- as.yearqtr(at)
    else x <- as.yearqtr(x)
    range <- par("usr")[if (side%%2) 
        1:2
    else 3:4]
    # range[1] <- ceiling(range[1])
    # range[2] <- floor(range[2])
    d <- range[2] - range[1]
    z <- c(range, x[is.finite(x)])
    class(z) <- "yearqtr"
    if (d > N1) { # axis has years only
        z <- structure(pretty(z), class = "yearqtr")
    } else if (d > N2) { # axis has all years and unlabelled quarters
        z <- seq(min(x), max(x), 0.25)
	# z <- seq(floor(min(x)), ceiling(max(x)))
    } else { # years and quarters
        z <- seq(min(x), max(x), 0.25)
    }
    if (!mat) 
        z <- x[is.finite(x)]
    z <- z[z >= range[1] & z <= range[2]]
    z <- sort(unique(z))
    class(z) <- "yearqtr"
    if (identical(labels, TRUE)) {
	if (missing(format)) format <- c("%Y", "Q%q")
	if (length(format) == 1) format <- c(format, "")
	if (d <= N2) labels <- format.yearqtr(z, format = format[2])
	idx <- format.yearqtr(z, format = "%q") == "1"
	labels <- rep(NA, length(z))
	labels[idx] <- format.yearqtr(z[idx], format = format[1])
    } else if (identical(labels, FALSE)) 
        labels <- rep("", length(z))
    axis(side, at = z, labels = labels, ...)
}

summary.yearqtr <- function(object, ...)
  summary(as.numeric(object), ...)

## convert from package date
as.yearqtr.date <- function(x, ...) {
	as.yearqtr(as.Date(x, ...))
}

mean.yearqtr <- function (x, ...) as.yearqtr(mean(unclass(x), ...))

Summary.yearqtr <- function (..., na.rm)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if (!ok) stop(.Generic, " not defined for yearqtr objects")
    val <- NextMethod(.Generic)
    class(val) <- oldClass(list(...)[[1]])
    val
}

Sys.yearqtr <- function() as.yearqtr(Sys.Date())

range.yearqtr <- function(..., na.rm = FALSE) {
	as.yearqtr(range.default(..., na.rm = na.rm))
}

unique.yearqtr <- function(x, incomparables = FALSE, ...) {
	as.yearqtr(unique.default(x, incomparables = incomparables, ...))
}

xtfrm.yearqtr <- function(x) as.numeric(x)
