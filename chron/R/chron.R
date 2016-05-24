"chron" <-
function(dates. = NULL, times. = NULL,
         format = c(dates = "m/d/y", times = "h:m:s"),
         out.format, origin.)
{
    if(is.null(format))
        format <- c(dates = "m/d/y", times = "h:m:s")
    if(missing(out.format)){
        if(is.character(format))
            out.format <- format
        else
            stop('must specify the "out.format" argument')
    }
    given <- c(dates = !missing(dates.), times = !missing(times.))
    if(is.null(default.origin <- getOption("chron.origin")))
        default.origin <- c(month = 1, day = 1, year = 1970)
    if(all(!given))
        ## dates and times missing
        return(structure(numeric(0),
                         format = format, origin = default.origin,
                         class = c("chron", "dates", "times")))
    if(inherits(dates., "dates")) {
        if(missing(origin.))
            origin. <- origin(dates.)
        else origin(dates.) <- origin.
    }
    else if(missing(origin.))
        origin. <- default.origin
    if(given["dates"] && !given["times"]) {
        ## presumably only dates
        if(missing(format) && inherits(dates., "dates"))
            format <- attr(dates., "format")
        fmt <- switch(mode(format),
                      character = ,
                      list = format[[1]],
                      name = ,
                      "function" = format,
                      NULL = c(dates = "m/d/y"),
                      stop("unrecognized format"))
        dts <- convert.dates(dates., format = fmt, origin. = origin.)
        tms <- dts - floor(dts)	
	## if dates include fractions of days create a full chron
        if(!all(is.na(tms)) && any(tms[!is.na(tms)] != 0))
            return(chron(dates. = floor(dts), times. = tms, format
                         = format, out.format = out.format, origin. = 
                         origin.))
        ofmt <- switch(mode(out.format),
                       character = ,
                       list = out.format[[1]],
                       name = ,
                       "function" = out.format,
                       NULL = c(dates = "m/d/y"),
                       stop("invalid output format"))
        attr(dts, "format") <- ofmt
        attr(dts, "origin") <- origin.
        class(dts) <- c("dates", "times")
        names(dts) <- names(dates.)
        return(dts)
    }
    if(given["times"] && !given["dates"]) {
        ## only times
        if(missing(format) && inherits(times., "times")) {
            format <- attr(times., "format")
            if(!is.name(format))
                format <- rev(format)[[1]]
        }
        fmt <- switch(mode(format),
                      character = ,
                      list = rev(format)[[1]],
                      name = ,
                      "function" = format,
                      NULL = c(times = "h:m:s"),
                      stop("invalid times input format"))
        tms <- convert.times(times., fmt)
        ofmt <- switch(mode(out.format),
                       character = ,
                       list = rev(out.format)[[1]],
                       name = ,
                       "function" = out.format,
                       NULL = c(dates = "m/d/y"),
                       stop("invalid times output format"))
        attr(tms, "format") <- ofmt
        class(tms) <- "times"
        names(tms) <- names(times.)
        return(tms)
    }
    ## both dates and times 
    if(length(times.) != length(dates.)) {
        if(length(times.) == 1)
            times. <- rep.int(times., length(dates.))
        else if(length(dates.) == 1)
            dates. <- rep.int(dates., length(times.))
        else
            stop(paste(deparse(substitute(dates.)), "and",
                       deparse(substitute(times.)), "must have equal lengths"))
    }
    if(missing(format)) {
        if(is.null(fmt.d <- attr(dates., "format")))
            fmt.d <- format[1]
        if(is.null(fmt.t <- attr(times., "format")))
            fmt.t <- format[2]
        if(mode(fmt.d) == "character" && mode(fmt.t) == "character")
            format <- structure(c(fmt.d, fmt.t),
                                names = c("dates", "times"))
        else {
            fmt.d <- if(is.name(fmt.d)) fmt.d else fmt.d[[1]]
            fmt.t <- if(is.name(fmt.t)) fmt.t else rev(fmt.t)[[1]]
            format <- list(dates = fmt.d, times = fmt.t)
        }
    }
    if(any(length(format) != 2, length(out.format) != 2))
        stop("misspecified chron format(s) length")
    if(all(mode(format) != c("character", "list")))
        stop("misspecified input format(s)")
    if(all(mode(out.format) != c("list", "character")))
        stop("misspecified output format(s)")
    dts <- convert.dates(dates., format = format[[1]], origin. = origin.)
    tms <- convert.times(times., format = format[[2]])
    x <- unclass(dts) + unclass(tms)
    attr(x, "format") <- out.format
    attr(x, "origin") <- origin.
    class(x) <- c("chron", "dates", "times")
    nms <- paste(names(dates.), names(times.))
    if(length(nms) && any(nms != ""))
        names(x) <- nms
    return(x)
}

as.chron <- function(x, ...) UseMethod("as.chron")
as.chron.default <- function (x, format, ...)
{
    if(inherits(x, "chron"))
        return(x)
    if(is.numeric(x)) {
        if (missing(format) || is.null(format)) return(chron(x, ...))
        else return(as.chron(as.POSIXct(format(x, scientific = FALSE),
                                        tz = "GMT", format = format),
                             ...))
	}
    if (is.character(x)) {
        if (missing(format) || is.null(format)) {
            out <- suppressWarnings(try(chron(x, ...), silent = TRUE))
            ## If this fails, try Date or datetime.
            if(inherits(out, "try-error")) {
                xx <- sub("T", " ", x)
                out <- if(!any(grepl(" ", x, fixed = TRUE)))
                    as.chron(as.Date(xx), ...)
                else
                    as.chron(as.POSIXct(xx, tz = "GMT"), ...)
            }
        } else {
            out <- as.chron(as.POSIXct(x, format = format, tz = "GMT"),
                            ...)
        }
        return(out)
    }
    stop("'x' cannot be coerced to a chron object")
}
as.chron.POSIXt <- function(x, offset = 0, tz = "GMT", ...)
{
    ## offset is in hours relative to GMT
    if(!inherits(x, "POSIXt")) stop("wrong method")
    x <- as.numeric(as.POSIXct(as.character(x, tz = tz), tz = "GMT")) +
        60 * round(60 * offset)
    tm <- x %% 86400
    # if(any(tm != 0))
        chron(dates. = x %/% 86400, times. = tm / 86400, ...)
    # else
    #   chron(dates. = x %/% 86400, ...)
}
as.chron.Date <- function(x, ...)
{
    chron(unclass(x), ...)
}

asChronYearFreq <-
function(x, frac = 0, holidays = FALSE, frequency, ...)
{
    stopifnot(isTRUE((12 / frequency) %% 1 == 0))
    x <- unclass(x)
    year <- floor(x + 0.001)
    month <- floor(12 * (x - year) + 1 + 0.5 + 0.001)
    dd.start <- as.Date(paste(year, month, 1, sep = "-"))
    nd <- 32 * 12 / frequency 
    dd.end <- dd.start + nd - as.numeric(format(dd.start + nd, "%d"))
    if(identical(holidays, FALSE))
        chron(((1 - frac) * as.numeric(dd.start) +
               frac * as.numeric(dd.end)),
              ...)
    else
        chron(sapply(seq_along(x), function(i) {
            s <- unclass(seq(dd.start[i], dd.end[i], by = "days"))
            h <- if(isTRUE(holidays)) is.holiday(s) else is.holiday(s, holidays)
            ss <- s[!is.weekend(s) & !h]
            quantile(ss, probs = frac, names = FALSE)
        }), ...)
}

as.chron.yearmon <-
function(x, frac = 0, holidays = FALSE, ...)
{
    asChronYearFreq(x, frac = frac, holidays = holidays,
                    frequency = 12, ...)
}

as.chron.yearqtr <-
function(x, frac = 0, holidays = FALSE, ...)
{
    asChronYearFreq(x, frac = frac, holidays = holidays,
                    frequency = 4, ...)
}

as.chron.ts <-
function(x, frac = 0, holidays = FALSE, ...)
{
    asChronYearFreq(time(x), frac = frac, holidays = holidays,
                    frequency = frequency(x), ...)
}

as.chron.factor <- function(x, ...) 
{
    as.chron(as.character(x), ...)
}

"is.chron" <-
function(x)
    inherits(x, "chron")

as.data.frame.chron <- as.data.frame.vector

"convert.chron" <-
function(x, format = c(dates = "m/d/y", times = "h:m:s"), origin.,
         sep = " ", enclose = c("(", ")"), ...)
{
    if(is.null(x) || !as.logical(length(x)))
        return(numeric(length = 0))
    if(is.numeric(x))
        return(x)
    if(!is.character(x) && all(!is.na(x)))
        stop(paste("objects", deparse(substitute(x)), 
                   "must be numeric or character"))
    if(length(format) != 2)
        stop("format must have length==2")
    if(missing(origin.)
       && is.null(origin. <- getOption("chron.origin")))
        origin. <- c(month = 1, day = 1, year = 1970)
    if(any(enclose != ""))
        x <- substring(x, first = 2, last = nchar(x) - 1)
    str <- unpaste(x, sep = sep)
    dts <- convert.dates(str[[1]], format = format[[1]],
                         origin. = origin., ...)
    tms <- convert.times(str[[2]], format = format[[2]], ...)
    dts + tms
}

"format.chron" <-
function(x, format = att$format, origin. = att$origin, sep = " ",
         simplify, enclosed = c("(", ")"), ...)
{
    att <- attributes(x)
    if(length(format) == 1L) {
        if(!nzchar(format))
            format <- "%Y-%m-%d %H:%M:%S"
        return(format(as.POSIXct(x), format = format, tz = "GMT"))
    }
    if(missing(simplify))
        if(is.null(simplify <- getOption("chron.simplify")))
            simplify <- FALSE
    dts <- format.dates(x, format[[1]], origin. = origin., simplify = 
                        simplify)
    tms <- format.times(x - floor(x), format[[2]], simplify = simplify)
    x <- paste(enclosed[1], dts, sep, tms, enclosed[2], sep = "")	
    ## output is a character object w.o class
    att$class <- att$format <- att$origin <- NULL
    attributes(x) <- att
    x
}

"new.chron" <-
function(x, new.origin = c(1, 1, 1970),
         shift = julian(new.origin[1], new.origin[2], new.origin[3],
         c(0, 0, 0)))
{
    cl <- class(x)
    class(x) <- NULL                    # get rid of "delim" attribute
    del <- attr(x, "delim")
    attr(x, "delim") <- NULL            # map formats
    format <- attr(x, "format")
    format[1] <- switch(format[1],
                        abb.usa = paste("m", "d", "y", sep = del[1]),
                        abb.world = paste("d", "m", "y", sep = del[1]),
                        abb.ansi = "ymd",
                        full.usa = "month day year",
                        full.world = "day month year",
                        full.ansi = "year month year",
                        format[1])
    if(length(format) == 2)
        format[2] <- switch(format[2],
                            military = "h:m:s",
                            format[2])
	attr(x, "format") <- format
    orig <- attr(x, "origin")
    if(is.null(orig)) {
        x <- x - shift
        attr(x, "origin") <- new.origin
    }
    ## (update origin after we assign the proper class!)
    ## deal with times as attributes 
    tms <- attr(x, "times")
    if(!is.null(tms)) {
        if(all(tms[!is.na(tms)] >= 1))
            tms <- tms/(24 * 3600)
        x <- x + tms
        class(x) <- c("chron", "dates", "times")
    }
    else class(x) <- c("dates", "times")
    x
}

print.chron <-
function(x, digits = NULL, quote = FALSE, prefix = "", sep = " ",
         enclosed = c("(", ")"), simplify, ...)
{
    if(!as.logical(length(x))) {
        cat("chron(0)\n")
        return(invisible(x))
    }
    if(missing(simplify) &&
       is.null(simplify <- getOption("chron.simplify")))
            simplify <- FALSE
    xo <- x
    x <- format.chron(x, sep = sep, enclosed = enclosed, simplify = 
                      simplify)
    print.default(x, quote = quote)
    invisible(xo)
}

unique.chron <-
function(x, incomparables = FALSE, ...) 
    x[!duplicated(x, incomparables, ...)]

xtfrm.chron <-
function(x)
    as.numeric(x)

pretty.chron <-
function(x, ...)
{
   if(!inherits(x, "times"))
       x <- chron(x)
   x <- as.POSIXct(x)
   attr(x, "tzone") <- "GMT"
   ans <- pretty(x, ...)
   structure(as.chron(ans), labels = attr(ans, "labels"))
}
