## class creation
yearmon <- function(x) structure(floor(12*x + .0001)/12, class = "yearmon")

## coercion to yearmon: always go via numeric
as.yearmon <- function(x, ...) UseMethod("as.yearmon")
as.yearmon.default <- function(x, ...) as.yearmon(as.numeric(x))
as.yearmon.numeric <- function(x, ...) yearmon(x)
as.yearmon.integer <- function(x, ...) structure(x, class = "yearmon")
as.yearmon.yearqtr <- function(x, frac = 0, ...) {
    if (frac == 0) yearmon(as.numeric(x)) else
    as.yearmon(as.Date(x, frac = frac), ...)
}
as.yearmon.dates <- 
as.yearmon.Date <- 
as.yearmon.POSIXt <- function(x, ...) as.yearmon(with(as.POSIXlt(x, tz="GMT"), 1900 + year + mon/12))
# as.jul.yearmon <- function(x, ...) jul(as.Date(x, ...)) # jul is from tis pkg
as.yearmon.mondate <-
as.yearmon.timeDate <-
as.yearmon.jul <- function(x, ...) as.yearmon(as.Date(x, ...))
as.yearmon.factor <- function(x, ...) as.yearmon(as.character(x), ...)
as.yearmon.character <- function(x, format = "", ...) {
   if (format == "") {
        nch <- nchar(gsub("[^-]", "", x))
		nch[is.na(x)] <- NA
		nch <- na.omit(nch)
        if (length(table(nch)) != 1) 
            stop("yearmon variable can only have one format")
		format <- if (all(nch == 0)) "%B %Y"
		else if (all(nch == 1)) "%Y-%m" else "%Y-%m-%d"
   }
   has.short.keys <- rep(regexpr("%[mbByY%]", format) > 0, length(x))
   has.no.others <- regexpr("%", gsub("%[mbByY%]", "", format)) < 0
   z <- ifelse(has.short.keys & has.no.others,
      as.Date( paste("01", x, sep = "-"), paste("%d", format, sep = "-"), ... ),
      as.Date(x, format, ...))
   as.yearmon(as.Date(z, origin = "1970-01-01"))
}
as.yearmon.ti <- function(x, ...) as.yearmon(as.Date(x), ...)

## coercion from yearmon
# returned Date is the fraction of the way through the period given by frac
as.Date.yearmon <- function(x, frac = 0, ...) {     
  x <- unclass(x)
  if(all(is.na(x))) return(as.Date(x))
  year <- floor(x + .001)
  ix <- !is.na(year)
  month <- floor(12 * (x - year) + 1 + .5 + .001)
  dd.start <- as.Date(rep(NA, length(year)))
  dd.start[ix] <- as.Date(paste(year[ix], month[ix], 1, sep = "-")) 
  dd.end <- dd.start + 32 - as.numeric(format(dd.start + 32, "%d"))
  as.Date((1-frac) * as.numeric(dd.start) + frac * as.numeric(dd.end), origin = "1970-01-01")
}
as.POSIXct.yearmon <- function(x, tz = "", ...) as.POSIXct(as.Date(x), tz = tz, ...)
as.POSIXlt.yearmon <- function(x, tz = "", ...) as.POSIXlt(as.Date(x), tz = tz, ...)
as.list.yearmon <- function(x, ...) lapply(seq_along(x), function(i) x[i])
as.numeric.yearmon <- function(x, ...) unclass(x)
as.character.yearmon <- function(x, ...) format.yearmon(x, ...)
as.data.frame.yearmon <- function(x, row.names = NULL, optional = FALSE, ...) 
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

## other methods for class yearmon
c.yearmon <- function(...)
    as.yearmon(do.call("c", lapply(list(...), as.numeric)))

cycle.yearmon <- function(x, ...) round(12 * (as.numeric(x) %% 1)) + 1

format.yearmon <- function(x, format = "%b %Y", ...) 
{
    if (length(x) == 0) return(character(0))
    xx <- format(as.Date(x), format = format, ...)
    names(xx) <- names(x)
    xx
}

print.yearmon <- function(x, ...) { 
    print(format(x), ...)
    invisible(x) 
}

months.yearmon <- function(x, abbreviate = FALSE) {
    months(as.Date(x), abbreviate = abbreviate)
}

quarters.yearmon <- function(x, abbreviate = FALSE) {
    quarters(as.Date(x), abbreviate = abbreviate)
}

"[.yearmon" <- function (x, ..., drop = TRUE) 
{
    cl <- oldClass(x)
    class(x) <- NULL
    val <- NextMethod("[")
    class(val) <- cl
    val
}

"[[.yearmon" <- function (x, ..., drop = TRUE) 
{
    cl <- oldClass(x)
    class(x) <- NULL
    val <- NextMethod("[[")
    class(val) <- cl
    val
}

MATCH.yearmon <- function(x, table, nomatch = NA, ...)
    match(floor(12*as.numeric(x) + .001), floor(12*as.numeric(table) + .001), nomatch = nomatch, ...)

Ops.yearmon <- function(e1, e2) {
    e1 <- as.numeric(as.yearmon(e1))
    e2 <- as.numeric(as.yearmon(e2))
    rval <- NextMethod(.Generic)
    if(is.numeric(rval)) rval <- yearmon(rval)
    return(rval)
}

"-.yearmon" <- function (e1, e2) 
{
    if (!inherits(e1, "yearmon")) 
        stop("Can only subtract from yearmon objects")
    if (nargs() == 1) 
	return(- as.numeric(e1))
    if (inherits(e2, "yearmon")) 
        return(as.numeric(e1) - as.numeric(e2))
    if (!is.null(attr(e2, "class"))) 
      stop("can only subtract yearmon objects and numbers from yearmon objects")
    yearmon(unclass(e1) - e2)
}

is.numeric.yearmon <- function(x) FALSE

Axis.yearmon <- function(x = NULL, at = NULL, ..., side, labels = NULL)
    axis.yearmon(x = x, at = at, ..., side = side, labels = TRUE)

axis.yearmon <- function (side, x, at, format, labels = TRUE, ..., N1 = 25, N2 = 2) {
    # If years in range > N1 then only years shown.  
    # If years in range > N2 then month ticks are not labelled.
    mat <- missing(at) || is.null(at)
    if (!mat) # at not missing
        x <- as.yearmon(at)
    else x <- as.yearmon(x)
    range <- par("usr")[if (side%%2) 
        1:2
    else 3:4]
    # range[1] <- ceiling(range[1])
    # range[2] <- floor(range[2])
    d <- range[2] - range[1]
    z <- c(range, x[is.finite(x)])
    class(z) <- "yearmon"
    if (d > N1) { # axis has years only
        z <- structure(pretty(z), class = "yearmon")
    } else if (d > N2) { # axis has all years and unlabelled months
        z <- seq(min(x), max(x), 1/12)
	# z <- seq(floor(min(x)), ceiling(max(x)))
    } else { # years and months
        z <- seq(min(x), max(x), 1/12)
    }
    if (!mat) 
        z <- x[is.finite(x)]
    z <- z[z >= range[1] & z <= range[2]]
    z <- sort(unique(z))
    class(z) <- "yearmon"
    if (identical(labels, TRUE)) {
	if (missing(format)) format <- c("%Y", "%b")
	if (length(format) == 1) format <- c(format, "")
	labels <- if (d <= N2) format.yearmon(z, format = format[2])
    else rep(NA, length(z))
	idx <- format.yearmon(z, format = "%m") == "01"
	labels[idx] <- format.yearmon(z[idx], format = format[1])
    } else if (identical(labels, FALSE)) 
        labels <- rep("", length(z))
    axis(side, at = z, labels = labels, ...)
}

summary.yearmon <- function(object, ...)
  summary(as.numeric(object), ...)

###

## convert from package date
as.yearmon.date <- function(x, ...) {
	as.yearmon(as.Date(x, ...))
}

mean.yearmon <- function (x, ...) as.yearmon(mean(unclass(x), ...))

Summary.yearmon <- function (..., na.rm)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if (!ok) stop(.Generic, " not defined for yearmon objects")
    val <- NextMethod(.Generic)
    class(val) <- oldClass(list(...)[[1]])
    val
}

Sys.yearmon <- function() as.yearmon(Sys.Date())

range.yearmon <- function(..., na.rm = FALSE) {
	as.yearmon(range.default(..., na.rm = na.rm))
}

unique.yearmon <- function(x, incomparables = FALSE, ...) {
	as.yearmon(unique.default(x, incomparables = incomparables, ...))
}

xtfrm.yearmon <- function(x) as.numeric(x)


