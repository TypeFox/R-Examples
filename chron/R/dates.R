## <FIXME>
## Need to improve consistency for preserving format and origin
## attributes ...
## </FIXME>

"dates"<-
function(x, ...)
{
    fmt <- attr(x, "format")
    ## Why?
    x <- chron(dates. = x, ...)
    ## <NOTE>
    ## This used to be floor.chron() ...
    ## </NOTE>
    cl <- oldClass(x)
    out <- floor(unclass(x))
    class(out) <- cl[!as.logical(match(cl, "chron", 0))]
    cl <- oldClass(x)
    attr(out, "format") <- fmt
    out
}

"Math.dates" <-
function(x, ...)
{
    ok <- switch(.Generic,
                 trunc = ,
                 round = ,
                 signif = ,
                 ceiling = ,
                 floor = TRUE,
                 FALSE)
    if(!ok)
        stop(paste(.Generic, "not defined for dates objects"))
    cl <- oldClass(x)
    class(x) <- NULL
    out <- NextMethod(.Generic)
    class(out) <- cl
    attr(out, "format") <- attr(x, "format")
    out
}

"Ops.dates" <-
function(e1, e2)
{
    ok <- switch(.Generic,
                 "+" = ,
                 "-" = ,
                 "<" = ,
                 ">" = ,
                 "==" = ,
                 "!=" = ,
                 "<=" = ,
                 ">=" = TRUE,
                 FALSE)
    if(nargs() == 1) {
        ## unary operators (only + is valid)
        if(.Generic == "+")
            return(e1)
        else
            stop(paste("unary", .Generic, "not defined for chron objects"))
    }
    if(!ok)
        stop(paste(.Generic, "not defined for chron objects"))
    dates.flg <- nchar(.Method)
    if(is.character(e1)) {
        e1 <- chron(e1, format = attr(e2, "format"), origin. = origin(e2))
        dates.flg[1] <- TRUE
    }
    if(is.character(e2)) {
        e2 <- chron(e2, format = attr(e1, "format"), origin. = origin(e1))
        dates.flg[2] <- TRUE
    }
    scalar <- !all(dates.flg)           # scalar operand?
    o1 <- origin(e1)
    o2 <- origin(e2)
    if(!scalar) {
        if(.Generic == "+")
            stop("chron objects may not be added together")
        if(any(o1 != o2)) {
            warning("different origin in dates arithmetic")
            origin(e2) <- o2 <- o1
        }
	}
    val <- NextMethod(.Generic)
    boolean <- match(.Generic, c("==", "!=", ">", ">=", "<", "<="),
                     nomatch = 0)
    if(boolean)
        return(val)                     # make sure origin wasn't dropped
    if(!inherits(val, "dates")) {
        attr(val, "origin") <- if(dates.flg[1]) o1 else o2
        class(val) <- unique(c(.Class, class(val)))
    }
    tms <- as.vector(val)
    tmp <- tms - floor(tms)	
    ## If a fractional scalar operand, then dates become chrons
    if(scalar && length(tmp <- tmp[!is.na(tmp)]) && any(tmp != 0)) {
        if(length(fmt.val <- attr(val, "format")) < 2)
            attr(val, "format") <- c(fmt.val, "h:m:s")
        class(val) <- c("chron", "dates", "times")
    }
    ## dates - dates is days
    if(!scalar && inherits(val, "dates")) {
        if(length(fmt.val <- attr(val, "format")) < 2)
            attr(val, "format") <- "h:m:s"
        else attr(val, "format") <- rev(attr(val, "format"))[[1]]
        attr(val, "origin") <- NULL
        val <- times(val)
    }
    val
}

"Summary.dates" <-
function(x, ...)
{
    ok <- switch(.Generic,
                 max = ,
                 min = ,
                 range = TRUE,
                 FALSE)
    if(!ok)
        stop(paste(.Generic, 
                   "not defined for objects that inherit from dates"))
    val <- NextMethod(.Generic)
    attr(val, "origin") <- origin(x)
    class(val) <- class(x)
    val
}

"[<-.dates" <-
function(x, ..., value)
{
    if(!as.logical(length(value)))
        return(x)                       # as per p.104 in the blue book
    if(!is.numeric(value) && !is.character(value) && !all(is.na(value)))
        stop("replacement of/with chron objects must be with times objects")
    ox <- origin(x)
    fmt <- attr(x, "format")
    if(!inherits(value, "dates"))
        value <- chron(value, format = fmt, origin. = ox)
    else if(any(ox != origin(value)))
        origin(value) <- ox
    cl <- oldClass(x)
    class(x) <- class(value) <- NULL
    x <- NextMethod(.Generic)
    attr(x, "format") <- fmt
    attr(x, "origin") <- ox
    class(x) <- cl
    x
}

"all.equal.dates" <-
function(..., tolerance = 1/(10 * 24 * 60 * 60))
    NextMethod("all.equal", ..., tolerance = tolerance)

as.data.frame.dates <- as.data.frame.vector

"c.dates" <-
function(..., recursive = FALSE)
{
    ## output will have the format and origin corresponding to the
    ## argument with earliest origin 
    dots <- list(...)
    is.dts <- unlist(lapply(dots, inherits, "dates"))
    o <- matrix(unlist(lapply(dots, origin)), nrow = 3)
    all.orig <- julian(o[1,  ], o[2,  ], o[3,  ], origin. = c(0, 0, 0))
    earliest <- min(all.orig)
    mdy <- month.day.year(earliest, origin. = c(0, 0, 0))
    orig <- c(mdy$month, mdy$day, mdy$year)
    n <- length(dots)
    fmt <- attr(dots[[(1:n)[is.dts][match(earliest, all.orig)]]], "format")
    out <- vector("list", length = n)
    for(i in 1:n) {
        x <- dots[[i]]	
	## note that NA's don't need any further processing
        if(!all(is.na(x))) {
            if(is.dts[i]) {
                if(any(origin(x) != orig))
                    origin(x) <- orig
            }
            else x <- chron(x, format = fmt, origin. = orig)
        }
        out[i] <- list(x)
    }
    out <- chron(unlist(out, use.names = FALSE),
                 origin. = orig, format = fmt)
    out
}

"convert.dates" <-
function(dates. = NULL, format = "m/d/y", origin., length. = 0, ...)
{
    ## returns a julian vector given various types of input
    if(is.null(dates.) || !length(dates.)) 
        return(numeric(length = length.))
    if(is.numeric(dates.))
        return(dates.)                  # assume julian format
    if(!is.character(dates.) && all(!is.na(dates.)))
        stop(paste("object", deparse(substitute(dates.)), 
                   "must be numeric or character"))
    if(!is.character(format)) {
        ## format may be a function or fun name
        FUN <- switch(mode(format),
                      name = get(as.character(format), mode = "function"),
                      "function" = format,
                      stop(paste("unrecognized date format",
                                 as.character(format))))
        return(FUN(dates., ...))
    }
    if(missing(origin.)
       && is.null(origin. <- getOption("chron.origin")))
        origin. <- c(month = 1, day = 1, year = 1970)	
    ## determine sep, order of month, day, year, etc.
    fmt <- parse.format(format)
    out <- if(nzchar(fmt$sep))
        unpaste(dates., sep = fmt$sep, fnames = fmt$periods,
                nfields = 3)
    else
        .str_to_ymd_list(dates., fmt)
    if(fmt$abb)
        mo <- as.numeric(out$m)
    else mo <- match(tolower(substring(out$m, 1, 3)),
                     tolower(month.abb), nomatch = NA)
    yy <- as.numeric(out$y)
    dy <- as.numeric(out$d)
    if(all(is.na(yy) | is.na(dy) | is.na(mo)))
        if(any(!is.na(as.character(dates.))))
            stop(paste("format", format, "may be incorrect"))
        else 
            return(rep(NA, length(dates.)))
    if(any(!is.na(yy)) && fmt$year.abb){
        fun <- getOption("chron.year.expand")
        fun <- switch(mode(fun), 
                      "character" = get(fun, mode = "function"),
                      "name" = eval(fun),
                      "function" = fun,
                      stop(paste("cannot expand 2-digit year abbreviation",
                                 "--you must specify \"chron.year.expand\"",
                                 "through options()")))
        yy <- fun(yy, ...)
    }
    non.na <- !is.na(mo)                # all months between 1 and 12?
    bad <- seq_along(mo)[non.na][mo[non.na] < 1 | mo[non.na] > 12]
    if(n.bad <- length(bad)) {
        if(n.bad > 10)
            msg <- paste(n.bad, "months out of range set to NA")
        else msg <- paste("month(s) out of range in positions",
                          paste(bad, collapse = ","), "set to NA")
        warning(msg)
        mo[bad] <- NA
        non.na[bad] <- FALSE
    }
    non.na <- non.na & !is.na(dy)
    mon.len <- month.length[mo[non.na]]
    mon.len[leap.year(yy[non.na]) & mo[non.na] == 2] <- 29# leap years!
    ## all days in the proper range (including leap years)?
    bad <- seq_along(dy)[non.na][dy[non.na] < 1 | dy[non.na] > mon.len]
    if(n.bad <- length(bad)) {
        if(n.bad > 10)
            msg <- paste(n.bad, "days out of range set to NA")
        else msg <- paste("days(s) out of range in positions", 
                          paste(bad, collapse = ","), "set to NA")
        warning(msg)
        dy[bad] <- NA
        non.na[bad] <- FALSE
    }
    return(julian(mo, dy, yy, origin. = origin.))
}

"cut.dates"<-
function(x, breaks, labels, start.on.monday = TRUE, ...)
{
    if(!inherits(x, "dates"))
        x <- chron(x)
    n <- length(breaks)                 # dates breaks may be either
                                        # numeric of character
    if(n > 1) {
        if(!inherits(breaks, "dates"))
            breaks <- sort(chron(dates. = breaks))	
	## make sure x and breaks have same origin
        org <- origin(x)
        if(!is.null(o <- origin(breaks)) && any(o != org))
            origin(breaks) <- org
        breaks <- as.numeric(breaks)
        if(missing(labels))
            labels <- paste("Range", seq_along(breaks[-1]))
        out <- cut.default(x, breaks = breaks, labels = labels)
        out <- ordered(as.character(out), levels = levels(out),
                       labels = labels)
        return(out)
    }
    if(n < 1) stop(paste(deparse(substitute(breaks)), 
                         "must have length > 0"))	
    ## breaks is either number or a string
    if(is.numeric(breaks)) {
        x <- as.numeric(x)
        if(inherits(breaks, "times"))
            breaks <- unclass(breaks)
        out <- NextMethod("cut")
        return(ordered(out))
    }
    ## we have a character string 
    valid <- c("days", "weeks", "months", "years")
    if(!as.logical(i <- pmatch(breaks[1], valid, 0)))
        stop(paste("unrecognized time period (", breaks, 
                   "), must be one of", paste(valid, collapse = ","), 
                   collapse = " "))
    by <- valid[i]
    bump <- c(1, 7, 31, 365)[i]         # force a full period for last obs.
    from <- min(x)
    orig <- origin(x)
    mdy <- month.day.year(as.numeric(from), origin. = orig)
    from <- switch(by,
                   days = from,
                   weeks = (from - day.of.week(mdy$m, mdy$d, mdy$y)
                            + as.numeric(start.on.monday)),
                   months = chron(julian(mdy$m, 1, mdy$y, origin. = orig)),
                   years = chron(julian(1, 1, mdy$y, origin. = orig)))
    if(from == min(x))
        from <- from - .Machine$double.eps
    breaks <- brk <- seq(from = from, to = max(x) + bump, by = by)
    breaks <- as.numeric(breaks)
    n <- length(breaks)
    x <- as.numeric(x)
    if(missing(labels)) {
        labels <-
            switch(by,
                   days = paste("day", seq_along(breaks[ - n] + 1)),
                   weeks = paste("week", seq_along(breaks[ - n] + 1)),
                   months = paste(as.character(months(brk[ - n] + 1)), 
                   substring(as.character(years(brk[ - n] + 1)), 3, 4)),
                   years = substring(as.character(years(brk[ - n] + 1)), 3, 4))
    }
    out <- cut.default(x, breaks = breaks, labels = labels, right = FALSE)
    ordered(as.character(out), levels = levels(out), labels = labels)
}

"format.dates" <-
function(x, format = "m/d/y", origin., simplify = FALSE, ...)
{
    if(!all(is.na(x)) && !is.numeric(x))
        stop(paste("couldn't extract julian dates from object", 
                   deparse(substitute(x))))
    if(is.null(default.orig <- getOption("chron.origin")))
        default.orig <- c(month = 1, day = 1, year = 1970)
    att <- attributes(x)
    if(inherits(x, "dates")) {
        if(missing(format))
            format <- switch(mode(att$format),
                             character = ,
                             list = att$format[[1]],
                             name = ,
                             "function" = att$format,
                             NULL = format,
                             stop("invalid output format for dates"))
        if(missing(origin.))
            origin. <- att$origin
    }
    else if(missing(origin.))
        origin. <- default.orig
    if(!is.character(format)) {
        ## format may be a function
        FUN <- switch(mode(format),
                      "function" = format,
                      name = eval(format),
                      stop(paste("unknown date format",
                                 as.character(format))))
        return(FUN(unclass(x), ...))
    }
    v <- month.day.year(floor(unclass(x)), origin. = origin.)
    v$day <- substring(paste("0", v$day, sep = ""), 
                       first = nchar(paste(v$day)))
    if(simplify) {
        drop.year <- length(unique(v$year[!is.na(v$year)])) <= 1
        drop.mon <- (simplify > 1 && drop.year
                     && length(unique(v$month)) <= 1)
        if(!drop.mon && !drop.year)
            drop.day <- TRUE
    }
    fmt <- parse.format(format[1])
    perm <- fmt$periods
    if(fmt$abb) {
        v$month <- substring(paste("0", v$month, sep = ""), 
                             first = nchar(paste(v$month)))
        if(fmt$year.abb){
            v$year <- v$year %% 100
            v$year <- substring(paste("0", v$year, sep=""),
                                first = nchar(paste(v$year)))
        }
    }
    else {
        v$month <- if(fmt$mon.abb)
            month.abb[v$month]
        else
            month.name[v$month]
    }
    sep <- fmt$sep
    y <- character(length = length(x))
    if(!simplify) {
        ## Perform partial matching by hand:
        ind <- pmatch(perm, names(v))
        y[] <- paste(v[[ind[1]]], v[[ind[2]]], v[[ind[3]]], sep = sep)
        ## "Simpler" than
        ##     do.call("paste",
        ##             c(v[pmatch(perm, names(v))], list(sep = sep))
        ## Could also use [[ with exact = FALSE, of course (R >= 2.6.0).
    } else {
        ## simplify (drop year/month when all equal)
        if(drop.mon) y[] <- v$day else if(drop.year) {
            perm <- perm[perm != "y"]	# drop years
            ind <- pmatch(perm, names(v))
            y[] <- paste(v[[ind[1]]], v[[ind[2]]], sep = sep)
        }
        else {
            perm <- perm[perm != "d"]	# drop days
            ind <- pmatch(perm, names(v))
            y[] <- paste(v[[ind[1]]], v[[ind[2]]], sep = sep)
        }
    }
    y[is.na(x)] <- NA
    y[x == Inf] <- "Inf"
    y[x ==  - Inf] <- "-Inf"
    att$format <- att$origin <- att$class <- NULL
    attributes(y) <- att
    y
}

print.dates <-
function(x, digits = NULL, quote = FALSE, prefix = "", simplify, ...)
{
    if(!as.logical(length(x))) {
        cat("dates(0)\n")
        return(invisible(x))
    }
    if(missing(simplify) &&
       is.null(simplify <- getOption("chron.simplify")))
            simplify <- FALSE
    print.default(format.dates(x, simplify = simplify), quote = quote)
    invisible(x)
}

seq.dates <- function(from, to, by = "days", length., ...)
{
    if(missing(from))
        stop("argument \"from\" must be specified")
    if(!inherits(from, "dates")) from <- chron(from[1])	
    ## the output will have same format and origin as "from"
    fmt <- attr(from, "format")         # dates format 
    org <- origin(from)                 # dates origin
    if(is.numeric(by)) {
        cl <- class(from)
        from <- as.numeric(from)
        if(!missing(to)) {
            if(!inherits(to, "dates")) to <- chron(to[1])
            if(!is.null(to.org <- origin(to)) && any(to.org != org))
                origin(to) <- org
            to <- as.numeric(to)
        }
        x <- seq.int(from, to, by)
	## preserve full chrons (i.e., don't round x)
        if(all(cl != "chron"))
            x <- round(x, 0)
        return(chron(x, format = fmt, origin. = org))
    }
    if(!is.character(by) || length(by) != 1)
        stop("\"by\" must be a number or string (days, weeks, months, or years)"
			)
    valid <- c("days", "weeks", "months", "years")
    if(!as.logical(i <- pmatch(by, valid, 0)))
        stop("\"by\" must be one of days, weeks, months, or years")
    by <- valid[i]                      # coerced "to" to a dates object
    if(missing(to)) {
        if(missing(length.))
            stop("must specify \"length\" when \"to\" is missing")
        to <- from + (length. - 1) * c(1, 7, 31, 366)[i]	
	## possibly BUGGY!!!
    }
    else {
        if(!missing(by) && !missing(length.))
            stop("Too many arguments")
        if(!inherits(to, "dates"))
            to <- chron(to)
        if(!missing(length.))
            by <- if(from < to) as.numeric(to - from)/(length. - 1) else 0
    }
    ## make sure "from" and "to" have the same origin
    if(!is.null(to.org <- origin(to)) && any(to.org != org))
        origin(to) <- org
    if(from > to)
        stop("\"from\" must be a date before \"to\"")
    frm <- as.numeric(from)
    t0 <- as.numeric(to)
    frm.mdy <- month.day.year(frm, origin. = org)	
    ## the idea is to generate all days between "form" and "to", subset
    ## out the dates we need, and finally chron them.
    x <- seq.int(from = frm, to = t0)
    if(by == "weeks") {
        mdy <- month.day.year(x, origin. = org)
        mdy.dow <- day.of.week(mdy$month, mdy$day, mdy$year)
        frm.dow <- day.of.week(frm.mdy$month, frm.mdy$day, frm.mdy$year)
        x <- x[mdy.dow == frm.dow]
    }
    else if(by == "months") {
        ## be careful when "from" is in the tail of the month!
        nxt.day <- month.day.year(as.numeric(from + 1))$month
        end.of.the.month <- frm.mdy$month != nxt.day
        if(end.of.the.month) x <- c(x, x[length(x)] + 1)
        mdy <- month.day.year(x, origin. = org)
        dys <- mdy$day
        if(frm.mdy$day <= 28)
            x <- x[dys == frm.mdy$day]
        else if(end.of.the.month)
            x <- x[dys == 1] - 1
        else {
            ## 29th or 30th of one of the 31-day months
            x1 <- x[dys == frm.mdy$day]	# all but Feb!
            x2 <- x[mdy$month == 3 & dys == 1] - 1 # Feb
            ## <NOTE>
            ## Of course, leap years can have Feb 29, in which case we
            ## get common entries in x1 and x2 ... hence, unique().
            x <- sort(unique(c(x1, x2)))
            ## </NOTE>
        }
        ## simple case
        if(!missing(length.)) x <- x[seq_len(length.)]
    }
    else if(by == "years") {
        ## be careful when "from" is Feb 29 of a leap year
        mdy <- month.day.year(x, org)
        if(leap.year(frm.mdy$year) && frm.mdy$day == 29)
            x <- x[mdy$day == 1 & mdy$month == 3] - 1
        else
            x <- x[mdy$day == frm.mdy$day & mdy$month == frm.mdy$month]
        if(!missing(length.)) x <- x[seq_len(length.)]
    }
    ## The original code had just
    ##   return(chron(x, format = fmt, origin = org))
    ## As pointed out by Sebastian Luque <sluque@mun.ca>, this causes
    ## trouble in case we have 00:00:00 time components, as in this case
    ## chron() returns a dates-only object.  Hence:
    if(inherits(from, "chron"))         # a full chron ...
        chron(floor(x), x - floor(x), format = fmt, origin. = org)
    else
        return(chron(x, format = fmt, origin. = org))
}

unique.dates <-
function(x, incomparables = FALSE, ...) 
    x[!duplicated(x, incomparables, ...)]

xtfrm.dates <-
function(x)
    as.numeric(x)

## chron 'dates' objects: only dates
## (no times here because caught by 'chron' method)
pretty.dates <-
function(x, ...)
{
   if(!inherits(x, "times"))
       x <- chron(x)
   x <- as.Date(x)
   ans <- pretty(x, ...)
   structure(as.chron(ans), labels = attr(ans, "labels"))
}
