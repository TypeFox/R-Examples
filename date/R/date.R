as.date <- function(x, order = "mdy", ...) {
    if (inherits(x, "date")) x
    else if (is.character(x)) {
	order.vec <-
            switch(order,
                   "ymd" = c(1, 2, 3),
                   "ydm" = c(1, 3, 2),
                   "mdy" = c(2, 3, 1),
                   "myd" = c(2, 1, 3),
                   "dym" = c(3, 1, 2),
                   "dmy" = c(3, 2, 1),
                   stop("Invalid value for 'order' option"))
	nn <- length(x)
	temp <- .C("char_date",
                   as.integer(nn),
                   as.integer(order.vec),
                   as.character(x),
                   month =integer(nn),
                   day = integer(nn),
                   year = integer(nn),
                   PACKAGE = "date")
	month <- ifelse(temp$month < 1 | temp$month > 12, NA, temp$month)
	day   <- ifelse(temp$day == 0, NA, temp$day)
	year  <- ifelse(temp$year == 0, NA, temp$year)
	temp <- mdy.date(month, day, year, ...)
    }
    else if (is.numeric(x)) {
	temp <- floor(x)
	attr(temp, "class") <- "date"
	}
    else stop("Cannot coerce to date format")
    temp
}

is.date <- function(x)
    inherits(x, "date")

Ops.date <- function(e1, e2) {
    ## Certain operation yield a date, others just give a number.  In
    ## order to make plotting functions work well, we end up allowing
    ## most all numeric operations.
    if (missing(e2))
        stop("Unary operations not meaningful for dates")
    if (.Generic == "&" || .Generic== "|")
	stop(paste("'", .Generic, "' not meaningful for dates",
                   sep = ""))
    class(e1) <- NULL
    class(e2) <- NULL
    if (.Generic == "-") {
	if (.Method[2] == "" ) {
            ## subtract a constant from a date 
            e1 <- as.integer(e1 - e2)
            class(e1) <- "date"
            e1
        }
	else if ((.Method[1] == "Ops.date" && .Method[2] == "Ops.date") ||
		 (.Method[1] == ""))
            e1 - e2
	else
            ## date - factor should fail
            stop("Invalid operation for dates")
	}
    else if (.Generic == "+") {
	if (.Method[1] == "" || .Method[2]=="")  {
            ## add constant to a date
            e1 <- as.integer(e1 + e2);
            class(e1) <- "date"
            e1
        }
	else e1 + e2
	}
    else get(.Generic)(e1, e2)
}

Math.date <- function(...)
    stop("Invalid operation on dates")

Summary.date <- function (..., na.rm = FALSE) {
    ok <- switch(.Generic, min = , max = , range = TRUE, FALSE)
    if (!ok)
        stop(paste(.Generic, "not defined for dates"))
    as.date(NextMethod(.Generic))
}

"[.date" <- function(x, ..., drop = TRUE) {
    cl <- class(x)
    class(x) <- NULL
    x <- NextMethod("[")
    class(x) <- cl
    x
}

"[[.date" <- function(x, ..., drop = TRUE) {
    cl <- class(x)
    class(x) <- NULL
    x <- NextMethod("[[")
    class(x) <- cl
    x
}

as.character.date <- function(x, ...) {
    fun <- options()$print.date
    if (is.null(fun))
        date.ddmmmyy(x)
    else
        get(fun)(x)
}

as.data.frame.date <- as.data.frame.vector

as.vector.date <- function(x, mode = "any") {
    if (mode == "any")
        as.vector(as.numeric(x), mode)
    else if (mode == "character" || mode == "logical" || mode == "list")
        as.vector(as.character(x), mode)
    else as.vector(unclass(x), mode)
}
    
is.na.date <- function(x) {
    NextMethod(.Generic)
}

plot.date <- function(x, y, ..., axes, xaxt, xlab, ylab,
                      xlim = range(x, na.rm = TRUE),
                      ylim = range(y, na.rm = TRUE))
{
    if(missing(xlab))
        xlab <- deparse(substitute(x))
    if(missing(ylab))
        ylab <- deparse(substitute(y))
    class(x) <- NULL                    # after deparse(substitute())
    if(!missing(axes) && !axes)         # argument axes works
        plot(x, y, ..., axes = axes, xlab = xlab, ylab = ylab,
            xlim = xlim, ylim = ylim)
    else if(!missing(xaxt))
        plot(x, y, ..., xaxt = xaxt, xlab = xlab, ylab = ylab,
            xlim = xlim, ylim = ylim)
    else {
        plot(x, y, ..., xaxt = "n", xlab = xlab, ylab = ylab,
            xlim = xlim, ylim = ylim)
        x <- c(x[!is.na(x)], xlim)      # draws axis completely when
                                        # using xlim
        xd <- date.mdy(x)
        ## get default for n from par("lab")
        temp <- pretty(x, n = par("lab")[1])
        delta <- temp[2] - temp[1]
        if(delta < 1)
            temp <- seq(min(x), max(x), 1)
        else if(delta > 182) {
            temp <- xd$year + (x - mdy.date(1, 1, xd$year))/365
            ## get default for n from par("lab")
            temp <- pretty(temp, n = par("lab")[1]) 
            temp <- mdy.date(1, 1, floor(temp)) + floor((temp %% 1) * 365)
        }
        axis(1, temp, as.character.date(temp), ...)
    }
}

print.date <- function(x, quote, prefix, ...) {
    if (missing(quote))
        quote <- FALSE
    invisible(print(as.character(x), quote = quote))
}

summary.date <- function(object, ...) {
    y <- as.character(range(object), ...)
    names(y) <- c("First ", "Last  ")
    y
}

Axis.date <-
function(x = NULL, at = NULL, xlim=range(x, na.rm=TRUE), ..., 
         side, labels = NULL)
{
    if(!is.null(x)) {
        x <- c(x[!is.na(x)], xlim)
        xd <- date.mdy(x)
        temp <- pretty(x, n = par("lab")[1L])
        delta <- temp[2L] - temp[1L]
        if(delta < 1)
            temp <- seq(min(x), max(x), 1)
        else if(delta > 182) {
            temp <- xd$year + (x - mdy.date(1, 1, xd$year)) / 365
            temp <- pretty(temp, n = par("lab")[1L])
            temp <- mdy.date(1, 1, floor(temp)) + floor((temp %% 1) * 365)
        }
        axis(side = side, at = temp, labels = as.character.date(temp), ...)
    } else {
        axis(side = side, at = at, labels = labels, ...)
    }
}

mdy.date <- function(month, day, year, nineteen = TRUE, fillday = FALSE,
                     fillmonth = FALSE) {
    ## Get the Julian date, but centered a la SAS, i.e., Jan 1 1960 is
    ## day 0.  Algorithm taken from Numerical Recipies.
    temp <- any((month != trunc(month)) |
                (day != trunc(day)) |
                (year != trunc(year)))
    if (!is.na(temp) && temp) {
	warning("Non integer input values were truncated in mdy.date")
	month <- trunc(month)
	day <- trunc(day)
	year <- trunc(year)
    }
    if (nineteen)
        year <- ifelse(year < 100, year + 1900, year)

    ## Force input vectors to be the same length, but in a way that
    ## gives an error if their lengths aren't multiples of each other.
    temp <- numeric(length(month + day + year))
    month <- month + temp
    day   <- day + temp
    year  <- year + temp

    if (fillmonth) {
	temp <- is.na(month)
	month[temp] <- 7
	day[temp] <- 1
	}
    if (fillday) day[is.na(day)] <- 15

    month[month < 1 | month > 12] <- NA
    day[day < 1] <- NA
    year[year == 0] <- NA               # there is no year 0
    year <- ifelse(year < 0, year + 1, year)
    tyear<- ifelse(month > 2, year, year - 1)
    tmon <- ifelse(month > 2, month + 1, month + 13)

    julian <-
        trunc(365.25 * tyear) + trunc(30.6001 * tmon) + day - 715940
    ## Check for Gregorian calendar changeover on Oct 15, 1582
    temp <- trunc(0.01 * tyear)
    save <- ifelse(julian >= -137774,
                   julian + 2 + trunc(.25 * temp) - temp,
                   julian)

    ## Check for invalid days (31 Feb, etc.) by calculating the Julian
    ## date of the first of the next month
    year <- ifelse(month == 12, year+1, year)
    month<- ifelse(month == 12, 1, month + 1)
    day <- 1
    tyear<- ifelse(month > 2, year, year - 1)
    tmon <- ifelse(month > 2, month + 1, month + 13)
    julian <-
        trunc(365.25 * tyear) + trunc(30.6001 * tmon) + day - 715940
    temp <- trunc(0.01 * tyear)
    save2<- ifelse(julian >= -137774,
                   julian + 2 + trunc(.25 * temp) - temp,
                   julian)

    temp <- as.integer(ifelse(save2 > save, save, NA))
    attr(temp, "class") <- "date"
    temp
}

date.mdy <- function(sdate, weekday = FALSE) {
    ##  Return the month, day, and year given a julian date
    attr(sdate, "class") <- NULL        # Stop any propogation of methods
    sdate <- floor(sdate + 2436935)     # From SAS to Num Recipies base
                                        # point 
    wday <- as.integer((sdate + 1) %% 7 +1)
    temp <- ((sdate - 1867216) -.25) / 36524.25
    sdate <- ifelse(sdate >= 2299161,
                    trunc(sdate+ 1 +temp - trunc(.25 * temp)),
                    sdate)
    jb <- sdate + 1524
    jc <- trunc(6680 + ((jb - 2439870) - 122.1) / 365.25)
    jd <- trunc(365.25 * jc)
    je <- trunc((jb - jd)/ 30.6001)
    day <- (jb - jd) - trunc(30.6001 * je)
    month <- as.integer(ifelse(je > 13, je - 13, je - 1))
    year  <- as.integer(ifelse(month > 2, jc - 4716, jc - 4715))
    year  <- as.integer(ifelse(year <= 0, year - 1, year))
    if (weekday)
        list(month = month, day = day, year = year, weekday = wday)
    else
        list(month = month, day = day, year = year)
}

date.ddmmmyy <- function(sdate) {
    temp <- date.mdy(sdate)
    tyr <- ifelse(floor(temp$year/100) == 19,
                  temp$year-1900, temp$year)
    month <- month.abb[temp$month]
    ifelse(is.na(sdate), as.character(NA),
           paste(temp$day, month, tyr, sep = ""))
}

date.mmddyy <- function(sdate, sep = "/") {
    temp <- date.mdy(sdate)
    tyr <- ifelse(floor(temp$year / 100) == 19,
                  temp$year - 1900, temp$year)
    ifelse(is.na(sdate), as.character(NA),
           paste(temp$month, temp$day, tyr, sep = sep))
}

date.mmddyyyy <- function(sdate, sep = "/") {
    temp <- date.mdy(sdate)
    ifelse(is.na(sdate), as.character(NA),
           paste(temp$month, temp$day, temp$year, sep = sep))
}
