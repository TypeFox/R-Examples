".Holidays" <-
    structure(.Data = c(8035, 8180, 8220, 8285, 8365, 8394),
              format = structure(.Data = "m/d/y", .Names = "dates"),
              origin = structure(.Data = c(1, 1, 1970),
              .Names = c("month", "day", "year")),
              class = c("dates", "times"),
              .Names = c("New Year's Day", "Memorial Day",
              "Independence Day", "Labor Day", "Thanksgiving",
              "Christmas"))
"day.abb" <-
    c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
"day.name" <-
    c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
      "Saturday")
"month.length"<-
    c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

"days"<-
function(x)
{
    if(!inherits(x, "dates"))
        x <- as.chron(x)
    d <- month.day.year(floor(as.numeric(x)), origin. = origin(x))$day
    ## use paste to avoid bug in ordered() as in beta release 8/92
    d <- ordered(paste(d), paste(1:31))
    d
}
"hours"<-
function(x)
{
    if(!inherits(x, "times"))
        x <- as.chron(x)        
    x <- as.numeric(x)
    sec <- round(24 * 3600 * abs(x - floor(x)))    
    hh <- sec %/% 3600
    hh
}
"minutes"<-
function(x)
{
    if(!inherits(x, "times"))
        x <- as.chron(x)        
    x <- as.numeric(x)
    sec <- round(24 * 3600 * abs(x - floor(x)))
    hh <- sec %/% 3600
    mm <- (sec - hh * 3600) %/% 60
    mm
}
"seconds"<-
function(x)
{
    if(!inherits(x, "times"))
        x <- as.chron(x)
    x <- as.numeric(x)
    sec <- round(24 * 3600 * abs(x - floor(x)))
    hh <- sec %/% 3600
    mm <- (sec - hh * 3600) %/% 60
    ss <- trunc(sec - hh * 3600 - 60 * mm)
    ss
}

"quarters.default"<-
function(x, abbreviate = TRUE)
{
    if(!inherits(x, "dates"))
        if((is.character(x) || is.numeric(x)))
            x <- chron(x)
        else return(NULL)
    v <- month.day.year(floor(as.numeric(x)))$month
    out <- (v - 1) %/% 3 + 1
    lbl <- if(abbreviate)
        c("1Q", "2Q", "3Q", "4Q")
    else
        c("I", "II", "III", "IV")
    out <- lbl[out]
    ordered(out, levels = lbl, labels = lbl)
}
"months.default"<-
function(x, abbreviate = TRUE)
{
    if(!inherits(x, "dates"))
        if((is.character(x) || is.numeric(x)))
            x <- chron(x)
        else return(NULL)
    out <- month.day.year(as.numeric(x), origin. = origin(x))$month
    lbl <- if(abbreviate) month.abb else month.name
    out <- lbl[out]
    ordered(out, levels = lbl, labels = lbl)
}

"weekdays.default" <-
function(x, abbreviate = TRUE)
{
    if(!inherits(x, "dates"))
        if((is.character(x) || is.numeric(x)))
            x <- chron(x)
        else stop("x must inherit from dates")
    v <- month.day.year(as.numeric(x), origin. = origin(x))
    out <- day.of.week(v$month, v$day, v$year) + 1
    lbl <- if(abbreviate) day.abb else day.name
    out <- lbl[out]
    ordered(out, levels = lbl, labels = lbl)
}

"years" <-
function(x)
{
    if(!inherits(x, "dates"))
        x <- as.chron(x)
    y <- month.day.year(as.numeric(x), origin. = origin(x))$year
    y <- ordered(y)
    y
}


"clock2frac" <-
function(str)
{
    h <- as.numeric(substring(str, 1, 2))
    m <- as.numeric(substring(str, 4, 5))
    w <- substring(str, 6, 7)
    if(any(h < 0, h > 12, m < 0, m > 59))
        stop("misspecified time")
    pm <- w == "pm" | w == "PM"
    h[pm] <- h[pm] + 12
    f <- (h * 3600 + m * 60)/(24 * 3600)
    f
}

"count.events" <-
function(x, by)
    table(cut(x, breaks = by))

"count.fields.str" <-
function(str, sep = "")
{
    n <- length(str)
    white.space <- missing(sep) || sep == ""
    .C("cnt_flds_str",
       strings = as.character(str),
       nstrings = as.integer(n),
       sep = as.character(sep),
       white.space = as.integer(white.space),
       counts = integer(n),
       PACKAGE = "chron")$counts
}

"day.of.week" <-
function(month, day, year)
{
    ix <- year + trunc((month - 14)/12)
    jx <- (trunc((13 * (month + 10 - (month + 10) %/% 13 * 12) - 1)/5)
           + day + 77 + (5 * (ix - (ix %/% 100) * 100)) %/% 4
           + ix %/% 400 - (ix %/% 100) * 2)
    jx %% 7
}

"format<-" <-
function(x, ..., value)
    UseMethod("format<-")

"frac2clock" <-
function(f)
{
    sec.per.day <- 24 * 3600
    secs <- f * sec.per.day
    h <- secs %/% 3600
    m <- round((secs - h * 3600)/60, 0)
    i <- h >= 13
    h[i] <- h[i] - 12
    pm <- rep("am", length(f))
    i <- f > 0.5
    pm[i] <- "pm"
    m <- paste(m)
    i <- nchar(m) == 1
    m[i] <- paste("0", m[i], sep = "")
    h <- paste(h)
    i <- nchar(h) == 1
    h[i] <- paste("0", h[i], sep = "")
    paste(h, ":", m, pm, sep = "")
}


"is.holiday" <-
function(x, holidays)
{
    if(!inherits(x, "dates")) x <- as.chron(x)
    if(missing(holidays)) {
        if(exists(".Holidays"))
            holidays <- .Holidays
        else holidays <- NULL
    } else if (length(holidays) == 0) holidays <- NULL
    if (is.null(holidays)) return(rep(FALSE, length(x)))
    orig.x <- origin(x)
    if(!is.null(orig.h <- origin(holidays)) && any(orig.x != orig.h))
        origin(holidays) <- orig.x
    out <- match(floor(x), floor(holidays), 0)
    as.logical(out)
}

"is.weekend" <-
function(x)
{
    if(!inherits(x, "dates")) x <- as.chron(x)
    v <- month.day.year(as.numeric(x), origin. = origin(x))
    out <- day.of.week(v$month, v$day, v$year) + 1	
    ## recall out is between 1 (Sunday) and 7 (Saturday)
    out == 1 | out == 7
}

"julian.default" <-
function(x, d, y, origin., ...)
{
    only.origin <- all(missing(x), missing(d), missing(y))
    if(only.origin) x <- d <- y <- NULL	# return days since origin
    if(missing(origin.) || is.null(origin.))
        if(is.null(origin. <- getOption("chron.origin")))
            origin. <- c(month = 1, day = 1, year = 1970)
    nms <- names(d)
    xdy <- cbind(x, d, y)
    m <- c(origin.[1], xdy[, "x"])      # prepend month of new origin
    d <- c(origin.[2], xdy[, "d"])      # prepend day of new origin
    y <- c(origin.[3], xdy[, "y"])      # prepend year of new origin
    ##
    ## code from julian date in the S book (p.269)
    ##
    y <- y + ifelse(m > 2, 0, -1)
    m <- m + ifelse(m > 2, -3, 9)
    c <- y %/% 100
    ya <- y - 100 * c
    out <- ((146097 * c) %/% 4 + (1461 * ya) %/% 4
            + (153 * m + 2) %/% 5 + d + 1721119)
    ## now subtract the new origin from all dates
    if(!only.origin) {
        if(all(origin. == 0))
            out <- out[-1]
        else
            out <- out[-1] - out[1]
        ## orig according to S algorithm
    }
    names(out) <- nms
    out
}

"julian2mine" <-
function(x)
{
    v <- month.day.year(x)
    d <- as.character(v$day)
    i <- nchar(d) == 1
    d[i] <- paste("0", d[i], sep = "")
    paste(d, month.abb[v$month], v$year, sep = "")
}

"leap.year" <-
function(y)
{
    if(inherits(y, "dates"))
        y <- month.day.year(as.numeric(y), origin. = origin(y))$year
    y %% 4 == 0 & (y %% 100 != 0 | y %% 400 == 0)
}

"mine2julian" <-
function(str)
{
    d <- substring(str, 1, 2)
    m <- substring(str, 3, 5)
    y <- substring(str, 6, 9)
    m <- match(m, month.abb, nomatch = NA)
    julian(m, as.numeric(d), as.numeric(y))
}

"month.day.year" <-
function(jul, origin.)
{
    if (!inherits(jul, "dates")) jul <- as.chron(jul)
    if(missing(origin.) || is.null(origin.))
        if(is.null(origin. <- getOption("chron.origin")))
            origin. <- c(month = 1, day = 1, year = 1970)
    if(all(origin. == 0)) shift <- 0 else shift <- julian(origin. = origin.)
    ## relative origin
    ## "absolute" origin
    j <- as.integer(floor(jul)) + as.integer(shift)
    j <- j - 1721119
    y <- (4 * j - 1) %/% 146097
    j <- 4 * j - 1 - 146097 * y
    d <- j %/% 4
    j <- (4 * d + 3) %/% 1461
    d <- 4 * d + 3 - 1461 * j
    d <- (d + 4) %/% 4
    m <- (5 * d - 3) %/% 153
    d <- 5 * d - 3 - 153 * m
    d <- (d + 5) %/% 5
    y <- 100 * y + j
    y <- y + ifelse(m < 10, 0, 1)
    m <- m + ifelse(m < 10, 3, -9)
    list(month = m, day = d, year = y)
}

"my.axis" <-
function(x, simplify = TRUE, ...)
{
    ## put date labels in one line plus time lables on second line
    px <- pretty(x)
    xx <- chron(px, format = attr(x, "format"), origin. = origin(x))
    lbls <- format(xx, enclose = c("", ""), sep = "\n", simplify = simplify)
    axis(1, at = px, labels = lbls, ...)
    invisible(list(at = px, labels = lbls))
}

"origin" <-
function(x)
    attr(x, "origin")
"origin<-" <-
function(x, value)
{
    if (length(value) != 3 || any(is.na(value)))
        stop("origin must be a month, day, year vector")
    if (value[1] < 1 || value[1] > 12)
        stop("month out of range in origin")
    n <- month.length[value[1]] +
        as.numeric(value[1] == 2 && leap.year(value[3]))
    if (value[2] < 1 || value[2] > n)
        stop("day out of range in origin")
    cl <- class(x)
    class(x) <- NULL
    jval <- julian(value[1], value[2], value[3], origin. = c(0, 0, 0))	
    ## adjust days for new origin (new.x + new.o == old.x + old.o)
    if (!is.null(ox <- attr(x, "origin")))
        x <- x - jval + julian(ox[1], ox[2], ox[3], origin. = c(0, 0, 0))
    new.origin <- unlist(month.day.year(jval, origin. = c(0, 0, 0)))
    attr(x, "origin") <-
        structure(new.origin, names = c("month", "day", "year"))
    class(x) <- cl
    x
}

"parse.format" <-
function(format, year.abb = getOption("chron.year.abb"), ...)
{
    ## determine order of month, day, year or hour, min, secs
    abb <- TRUE                         # short notation?
    mon.abb <- FALSE                    # should month names be abbreviated?
    if(is.null(year.abb))
        year.abb <- TRUE
    if((nf <- nchar(format)) == 5) {
        ## abbreviated dates/times
        sep <- substring(format, 2, 2)
        fmt <- substring(format, first = c(1, 3, 5), last = c(1, 3, 5))
    }
    else if(nf == 3) {
        sep <- ""                       # no sep
        fmt <- substring(format, first = 1:3, last = 1:3)
    }
    else {
        ## full format (month names)
        abb <- FALSE
        sep <- gsub("^[[:alpha:]]+([^[:alpha:]]).*", "\\1", format)
        if(sep == format)
            stop(paste("unrecognized format", format))
        fmt <- unlist(unpaste(format, sep = sep))
        mon.abb <- if(any(fmt == "month")) FALSE else TRUE
    }
    periods <- substring(tolower(fmt), 1, 1) # m, d, & y in right order
    return(list(abb = abb, sep = sep, periods = periods, 
                mon.abb = mon.abb, year.abb = year.abb))
}

"unpaste" <-
function(str, sep = "/", fnames = NULL, nfields = NULL,
         first = c(1, 3, 5), width = 2)
{
    ## split str into fields separated by sep or by fiels specified by
    ## start positions and field widths; output a list 
    str <- as.character(str)
    nas <- is.na(str) | str == ""
    if(sep != "") {
        if(is.null(nfields)) {
            ## use a simple heuristic
            nf <- count.fields.str(str[!nas], sep = sep)
            cnt <- table(nf)
            nfields <- sort(unique(nf))[cnt == max(cnt)]
        }
        str[nas] <- paste(rep(NA, nfields), collapse = sep)
        nf <- count.fields.str(str, sep = sep)
        bad <- seq_along(str)[nf != nfields]
        if(n.bad <- length(bad)) {
            if(n.bad > 10)
                msg <- paste(n.bad, 
                             "entries set to NA",
                             "due to wrong number of fields")
            else msg <- paste(
                              "wrong number of fields in entry(ies)",
                              paste(bad, collapse = ", "))
            warning(msg)
            nas[bad] <- TRUE
            str[nas] <- paste(rep(NA, nfields), collapse = sep)
        }
        n <- length(str)
        white.space <- FALSE
        out <- .Call("unpaste",
                      as.character(str),
                      as.character(sep),
                      as.logical(white.space),
                      as.integer(nfields),
                      PACKAGE = "chron")
        for(i in seq_along(out))
            out[[i]][nas] <- as.character(NA)
    }
    else {
        last <- first + width - 1
        out <- vector("list", length = length(first))
        for(i in seq_along(first)) {
            out[[i]] <- substring(str, first[i], last[i])
            out[[i]][nas] <- as.character(NA)
        }
    }
    names(out) <- fnames
    return(out)
}

.str_to_ymd_list <-
function(str, fmt)
{
    str <- as.character(str)
    nas <- is.na(str) | str == ""
    periods <- fmt$periods
    widths <- cbind(y = nchar(str) - 4, m = 2, d = 2)
    last <- apply(widths[, fmt$periods, drop = FALSE], 1, cumsum)
    first <- rbind(0, last[-3, , drop = FALSE]) + 1
    out <- vector("list", length = 3)
    for(i in seq_along(periods)) {
        out[[i]] <- substring(str, first[i, ], last[i, ])
        out[[i]][nas] <- as.character(NA)
    }
    names(out) <- periods
    out
}
