## Align a dates on a day, bizday, month, week or year boundary.

dateAlign <- function(x, by = 'days', k.by = 1, direction = 1, week.align = NULL, holidays = NULL, silent = FALSE)
    UseMethod("dateAlign")

dateAlign.character <- function(x, by = 'days', k.by = 1, direction = 1, week.align = NULL, holidays = NULL, silent = FALSE)
{
    x <- NextMethod('dateAlign')
    as.character(x)
}

dateAlign.POSIXct <- function(x, by = 'days', k.by = 1, direction = 1, week.align = NULL, holidays = NULL, silent = FALSE)
{
    tz <- attr(date, 'tzone')
    x <- NextMethod('dateAlign')
    # need to convert Date to character before converting back to POSIXct
    # see examples in tests/gotchas.Rt
    x <- as.POSIXct(as.character(x))
    if (!is.null(tz))
        attr(x, 'tzone') <- tz
    return(x)
}

dateAlign.POSIXlt <- function(x, by = 'days', k.by = 1, direction = 1, week.align = NULL, holidays = NULL, silent = FALSE)
{
    tz <- attr(date, 'tzone')
    x <- NextMethod('dateAlign')
    # need to convert Date to character before converting back to POSIXlt
    # see examples in tests/gotchas.Rt
    x <- as.POSIXlt(as.character(x))
    if (!is.null(tz))
        attr(x, 'tzone') <- tz
    return(x)
}


dateAlign.Date <- function(x, by = 'days', k.by = 1, direction = 1, week.align = NULL, holidays = NULL, silent = FALSE)
{
    ### BEGIN ARGUMENT PROCESSING ###
    if (!inherits(x, "Date"))
    {
        x <- dateParse(x)
        if (is.null(x))
          stop("'x' argument must inherit from the 'Date' class.")
    }

    if (!is.character(by))
        stop("'by' must be a character vector.")

    if (length(by) > 1)
        stop("'by' must be scalar.")

    if (length(byhol <- strsplit(by, '@')[[1]]) == 2) {
        if (!is.null(holidays))
            stop("cannot have both holidays = ", holidays, " and by = '", by, "'")
        by <- byhol[1]
        holidays <- byhol[2]
    }

    if (!(by %in% c('days', 'bizdays', 'weeks', 'months', 'years')))
        stop("'by' must be one of 'days', 'bizdays', 'weeks', 'months', 'years'.")

    if (length(k.by) > 1)
        stop("'k.by' must be scalar.")

    if (!is.null(week.align) && by != 'weeks')
        warning("ignoring week.align = ", week.align, " when by != 'weeks'")

    k.by <- as.integer(k.by)

    if (by == "days" && !(k.by >= 1 && k.by <= 30))
        stop("when using 'by = \"days\"', 'k.by' must be in the range: 1:30.")
    else
        if (k.by < 1) stop("'k.by' must be positive.")

    ## 'direction' is +1 for after, -1 for before
    direction <- as.integer(direction)
    if (!(direction == -1 || direction == 1))
        stop("'direction' must be -1 or 1.")

    if (!is.null(week.align) && (!is.numeric(week.align) || (week.align < 0 || week.align > 6)))
        stop("'week.align' must be between 0 and 6, where 0 is Sunday.")

    if (!is.null(holidays))
    {
        if (by != 'bizdays')
        {
            if (!silent)
                warning("ignoring 'holidays' argument. Only relevant when 'by = \"bizdays\"'.")
            holidays <- NULL
        }
    }

    switch(by,
           bizdays = k.by <- 1,
           weeks = {
               k.by <- 1
               if(is.null(week.align))
                   by <- "days"
           }
           )

    ### END ARGUMENT PROCESSING ###

    ## ALIGNMENT GUIDE
    ##
    ## Alignment of dates can be thought of as a partition on date sequences
    ## where an input date is aligned to the first date in a partition, if it is not
    ## already aligned. The direction of alignment determines which partition to use for
    ## the alignment. If the direction is <0 then alignment happens in the partition which
    ## the date falls in. If >0 then alignment happens in the partition just after the
    ## partition in which the dates falls.
    ##
    ## In the following examples, the pipe character delimits partitions, the 'i'
    ## denotes the "input" dates, and the '*' denotes the aligned date or "output".
    ##
    ## WEEK ALIGNMENT
    ##
    ## 1. by='weeks',  week.align=2, meaning align to Tuesday, and direction=-1.
    ## The input date is '2007/12/06', a Thursday, and the aligned date is '2007/12/04':
    ##
    ##      *   i
    ## 2 3 |4 5 6 7 8 9 10 |11 12 13 14
    ##
    ## 2.  by='weeks',  week.align=2, meaning align to Tuesday, and direction=1.
    ## The input date is '2007/12/06', a Thursday, and the aligned date is '2007/12/11':
    ##
    ##          i           *
    ## 2 3 |4 5 6 7 8 9 10 |11 12 13 14
    ##
    ## 3.  by='weeks',  week.align=2, meaning align to Tuesday, and direction=-1|1.
    ## The input date is '2007/12/04', a Tuesday, and the aligned date is '2007/12/04':
    ##
    ##      i
    ##      *
    ## 2 3 |4 5 6 7 8 9 10 |11 12 13 14
    ## DAYS ALIGNMENT
    ##
    ## Partitions for alignment by days, when k.by>1, starts on the first day of the month of the input
    ## date. It's unclear if there's a real-world use for this, though, so anyone having a better
    ## idea of where to start the partition, please let me know. other ideas are the begining of
    ## the year, or the beginning of the week of the input date)
    ##
    ## 3. by='days',  k.by=3, and direction=1.
    ## The input date is '2007/12/06', a Thursday, and the aligned date is '2007/12/07':
    ##
    ##             i  *
    ## |1 2 3 |4 5 6 |7 8 9 |10 11 12 ...
    ##
    ## 4. by='days', k.by=30, and direction=1
    ## The input date is '2007/12/06', a Thursday, and the aligned date is '2008/01/01':
    ##

    len <- length(x)

    if (by == 'days' || by == 'bizdays')
    {
        x <- as.POSIXlt(x)

        align <- function(x)
        {
            day <- x$mday - 1

            if (day %% k.by != 0)
            {
                part <- day %/% k.by

                if (direction > 0)
                    part <- part + 1

                x$mday <- part * k.by + 1

                ## If we go beyond end of month start with first
                ## day of next month. As a test for going beyond
                ## end of month we pretend that months are 31
                ## days long (extecpt for Feb. where we use 29).
                ## This is off by one day for about half the months
                ## but that is OK. It is only a problem if we are
                ## off by more that one day!
                if (direction > 0 && ((x$mday > 31 || x$mon == 1 && x$mday > 29)))
                {
                    x$mday <- 1
                    x$mon <- x$mon + 1
                }
            }
            x
        }

        x <- do.call("combine", lapply(seq(1, len), function(i) align(x[i])))

        if (by == 'bizdays')
        {
            align <- function(x)
            {
                ## Move over weekends and holidays.
                wday <- x$wday

                while (wday == 6 || wday == 0 ||
                       (!is.null(holidays) && isHoliday(x, holidays)))
                {
                    ## Step one day.
                    x$mday <- x$mday + direction
                    wday <- (wday + direction) %% 7
                }
                x
            }

            x <- do.call("combine", lapply(seq(1, len), function(i) align(x[i])))
        }
    }
    else if (by == 'weeks' && !is.null(week.align))
    {
        xp <- as.POSIXlt(x)

        align2 <- function(x, xp)
        {
            ## We are done if weekday is already 'week.align'.
            if (week.align != xp$wday)
            {
                ## An example partition by week with the first day on 2 (tuesday)
                ## 0 1 |2 3 4 5 6 0 1 |2 3 4 5 6 0
                ## 'forward' is how many dates we need to move forward
                forward <- ((week.align - xp$wday - 1) %% 7) + 1
                x <- as.Date(x) + forward
                if (direction < 0 && forward > 0)
                    x <- x - 7
            }
            x
        }

        x <- do.call("c", lapply(seq(1, len), function(i) align2(x[i], xp[i])))
    }
    else if (by == 'months')
    {
        x <- as.POSIXlt(x)

        align <- function(x)
        {
            ## We are done if month is already dividable with
            ## 'k.by' and it is the first of the month.
            if (!(x$mon %% k.by == 0 && x$mday == 1))
            {
                x$mday <- 1

                part <- x$mon %/% k.by

                if (direction > 0)
                    part <- part + 1

                x$mon <- part * k.by
            }
            x
        }

        x <- do.call("combine", lapply(seq(1, len), function(i) align(x[i])))
    }
    else if (by == 'years')
    {
        x <- as.POSIXlt(x)

        align <- function(x)
        {
            ## Start partitioning at year zero. This may seem arbitrary, but that
            ## is what S-PLUS does. POSIXlt has origin at 1900/1/1.
            year <- x$year + 1900

            ## We are done if year is already dividable with
            ## 'k.by' and the date is '1/1'.
            if (!(year %% k.by == 0 && x$mday == 1 && x$mon == 0))
            {
                x$mday <- 1
                x$mon <- 0

                part <- year %/% k.by

                if (direction > 0)
                    part <- part + 1

                x$year <- part * k.by - 1900
            }
            x
        }

        x <- do.call("combine", lapply(seq(1, len), function(i) align(x[i])))
    }

    as.Date(x)
}

dateAlign.default <- dateAlign.Date
