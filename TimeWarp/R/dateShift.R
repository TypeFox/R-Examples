## Shift dates a number of days, bizdays, months, weeks or years.

dateShift <- function(x, by = 'days', k.by = 1, direction = 1, holidays = NULL, silent = FALSE)
    UseMethod("dateShift")

dateShift.character <- function(x, by = 'days', k.by = 1, direction = 1, holidays = NULL, silent = FALSE)
{
    x <- NextMethod('dateShift')
    as.character(x)
}

dateShift.POSIXct <- function(x, by = 'days', k.by = 1, direction = 1, holidays = NULL, silent = FALSE)
{
    tz <- attr(date, 'tzone')
    x <- NextMethod('dateShift')
    # need to convert Date to character before converting back to POSIXct
    # see examples in tests/gotchas.Rt
    x <- as.POSIXct(as.character(x))
    if (!is.null(tz))
        attr(x, 'tzone') <- tz
    return(x)
}

dateShift.POSIXlt <- function(x, by = 'days', k.by = 1, direction = 1, holidays = NULL, silent = FALSE)
{
    tz <- attr(date, 'tzone')
    x <- NextMethod('dateShift')
    # need to convert Date to character before converting back to POSIXlt
    # see examples in tests/gotchas.Rt
    x <- as.POSIXlt(as.character(x))
    if (!is.null(tz))
        attr(x, 'tzone') <- tz
    return(x)
}

dateShift.Date <- function(x, by = 'days', k.by = 1, direction = 1, holidays = NULL, silent = FALSE)
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

    k.by <- as.integer(k.by)

    if (by == "days" && !(k.by >= 1 && k.by <= 30))
        stop("when using 'by = \"days\"', 'k.by' must be in the range: 1:30.")
    else
        if (k.by < 1) stop("'k.by' must be positive.")

    ## 'direction' is +1 for after, -1 for before
    direction <- as.integer(direction)
    if (!(direction == -1 || direction == 1))
        stop("'direction' must be -1 or 1.")

    if (!is.null(holidays))
    {
        if (by != 'bizdays')
        {
            if (!silent)
                warning("ignoring 'holidays' argument. Only relevant when 'by = \"bizdays\"'.")
            holidays <- NULL
        }
    }

    ### END ARGUMENT PROCESSING ###

    len <- length(x)
    x <- as.POSIXlt(x)

    if (by == 'days')
    {
        x$mday <- x$mday + k.by * direction
    }
    else if (by == 'bizdays')
    {
        shift <- function(x)
        {
            move <- 0

            ## Keep track of current weekday.
            wday <- x$wday

            while (move < k.by)
            {
                ## Step one day.
                x$mday <- x$mday + direction
                wday <- (wday + direction) %% 7

                ## Do not count weekend days and
                ## holidays.
                if (wday != 6 && wday != 0 &&
                    (is.null(holidays) || !isHoliday(x, holidays)))
                    move <- move + 1
            }
            x
        }

        x <- do.call("combine", lapply(seq(1, len), function(i) shift(x[i])))
    }
    else if (by == 'weeks')
    {
        x$mday <- x$mday + 7 * k.by * direction
    }
    else if (by == 'months')
    {
        x$mon <- x$mon + k.by * direction
    }
    else if (by == 'years')
    {
        x$year <- x$year + k.by * direction
    }

    as.Date(x)
}

dateShift.default <- dateShift.Date
