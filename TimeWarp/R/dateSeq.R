dateSeq <- function(from = NULL, to = NULL, year = NULL, by = "days",
                    k.by = 1, length.out = NULL, holidays = NULL,
                    align.by = TRUE, extend = FALSE, range = NULL,
                    week.align = NULL)
    UseMethod("dateSeq")

dateSeq.character <- function(from = NULL, to = NULL, year = NULL, by = "days",
                    k.by = 1, length.out = NULL, holidays = NULL,
                    align.by = TRUE, extend = FALSE, range = NULL,
                    week.align = NULL)
{
    x <- NextMethod('dateSeq')
    as.character(x)
}

dateSeq.POSIXct <- function(from = NULL, to = NULL, year = NULL, by = "days",
                    k.by = 1, length.out = NULL, holidays = NULL,
                    align.by = TRUE, extend = FALSE, range = NULL,
                    week.align = NULL)
{
    tz <- attr(date, 'tzone')
    x <- NextMethod('dateSeq')
    # need to convert Date to character before converting back to POSIXct
    # see examples in tests/gotchas.Rt
    x <- as.POSIXct(as.character(x))
    if (!is.null(tz))
        attr(x, 'tzone') <- tz
    return(x)
}

dateSeq.POSIXlt <- function(from = NULL, to = NULL, year = NULL, by = "days",
                    k.by = 1, length.out = NULL, holidays = NULL,
                    align.by = TRUE, extend = FALSE, range = NULL,
                    week.align = NULL)
{
    tz <- attr(date, 'tzone')
    x <- NextMethod('dateSeq')
    # need to convert Date to character before converting back to POSIXlt
    # see examples in tests/gotchas.Rt
    x <- as.POSIXlt(as.character(x))
    if (!is.null(tz))
        attr(x, 'tzone') <- tz
    return(x)
}


dateSeq.Date <- function(from = NULL, to = NULL, year = NULL, by = "days",
                    k.by = 1, length.out = NULL, holidays = NULL,
                    align.by = TRUE, extend = FALSE, range = NULL,
                    week.align = NULL) {
    ## Checking 'by' and 'holidays' arguments.
    if ((at.pos <- regexpr("@",by)[1]) != -1) {
        if (length(holidays)!=0)
            stop("cannot supply both holidays= and by=\"x@holidays\"")
        holidays <- substring(by, at.pos+1)
        if (holidays=="")
            stop("could not parse holiday name out of '", by, "'")
        by <- substring(by, 1, at.pos-1)
    }

    if ((space.pos <- regexpr(" ",by,)[1]) != -1) {
        a.str <- substring(by, 1, space.pos-1)
        by <- substring(by, space.pos+1)
        k.by <- as.numeric(a.str)
        if (is.na(k.by)) stop("k.by must be numeric")
    }

    ## Do we need these? I use the first one to match on 'by' below.
    by.choices.cont <- c('days', 'bizdays','weeks', 'months', 'years')
    by.choices.disc <- c('quarters','weekdays')
    if(!(by %in% by.choices.cont))
        stop(paste("by must be one of:",paste(by.choices.cont,collapse=","),"."))

    ## return a sequence of business days
    if (!is.null(range)) {
        if (!is.null(from))
            stop("cannot supply 'from' and 'range'")
        if (!is.null(to))
            stop("cannot supply 'to' and 'range'")
        if (!inherits(range, "Date"))
            range <- range(dateParse(range))
        from <- min(range)
        to <- max(range)
    }

    if (!is.null(year) && !is.null(from) && !is.null(to))
        stop("must supply either 'year' or 'from'/'to'")
    if (!is.null(length.out) && sum(is.null(from),is.null(to),is.null(year))!=2)
        stop("must supply one of 'from'/'to'/'year' with 'length.out'")
    if (is.null(year) && sum(is.null(from),is.null(to),is.null(length.out))!=1)
        stop("must supply 'year' or two of 'from', 'to' and 'length.out'")
    if (!is.null(year)) {
        if (is.null(from))
            from <- dateParse(paste(min(year), "/01/01", sep=''), format="%Y/%m/%d")
        if (is.null(to) && is.null(length.out))
            to <- dateParse(paste(max(year), "/12/31", sep=''),format="%Y/%m/%d")
    } else {
        if (!is.null(from) && !inherits(from, "Date"))
            from <- dateParse(from)
        if (!is.null(to) && !inherits(to, "Date"))
            to <- dateParse(to)
    }

    if (is.character(holidays)) {
        ## from and to are non-null if either they were specified or year was specified
        ## can only be null if they were not specified and year was non specified, which
        ## means that length.out was specified, and the other of from/to was specified
        if (by=="years")
            holidays.pad <- length.out+2
        else
            holidays.pad <- 2*length.out+10

        if (!is.null(from))
            from.bound <- dateWarp(from,-2,by=by)
        else
            from.bound <- dateWarp(to,-holidays.pad,by=by)
        if (!is.null(to))
            to.bound <- dateWarp(to,2,by=by)
        else
            to.bound <- dateWarp(from,holidays.pad,by=by)
        holidays <- holidays(years(from.bound) : years(to.bound), holidays)
    }

    ## Alignment
    if (!is.null(week.align) && (!is.numeric(week.align) || (week.align < 0 || week.align > 6)))
        stop("week.align must be between 0 and 6, where 0 is Sunday")

    if (!is.null(to) && !is.null(from) && from > to)
    {
        warning("'from' date is later in time than the 'to' date.")
        return(emptyDate())
    }

    align.dir <- if (extend) -1 else 1
    if (align.by)
    {
        ## Don't use holidays for align.by
        if (!is.null(from))
            from <- dateAlign(from, by=by, k.by=k.by, week.align=week.align, direction=align.dir)

        if (!is.null(to))
            to <- dateAlign(to, by=by, k.by=k.by, week.align=week.align, direction= -align.dir)
    }

    if (!is.null(to) && !is.null(length.out)) {
        from1 <- dateShift(to, by=by, k.by=k.by * (length.out-1), direction=-1, holidays=holidays)
        if (align.by)
            from1 <- dateAlign(from1, by=by, k.by=k.by, week.align=week.align, direction=align.dir)
        if (!is.null(from) && from1 != from)
            stop("supplied 'from'=", from, " is inconsistent with ", from1, " calculated from 'to' and 'length.out'")
        from <- from1
    }
    if ((!is.null(from)) + (!is.null(to)) + (!is.null(length.out)) < 2)
        stop("must specify two of 'from', 'to' and 'length.out'")

    ## Might have made from > to by the alignment -- don't warn about this
    if (!is.null(to) && from > to)
        return(emptyDate())

    ### END of argument checking

    x <- NULL
    if (by != 'bizdays'){
        rby <- paste(k.by,sub('s','',by))

        ## Don't supply 'to' because we always calculate 'from'
        if (!is.null(from) && !is.null(to))
            x <- seq(from=from, to=to, by=rby)
        else
            x <- seq(from=from, by=rby,length.out=length.out)
    } else if (by == 'bizdays'){
        if (is.null(to)){
            x <- rep(from,length.out)
            i <- 1
            w <- as.POSIXlt(from)$wday
            while(length.out){
                if ( (is.null(holidays) || !any(from==holidays)) && w != 0 && w != 6){
                    ## Day is neither a weekend day nor is it a holiday
                    ## so count as a business day
                    x[i] <- from
                    i <- i + 1
                    length.out <- length.out - 1
                }
                ## Increment/Decrement date to next day
                from <- from + 1

                ## Increment to next day, wrapping by 7
                w <- (w + 1) %% 7
            }
        } else {
            if (!is.null(length.out)) x <- rep(from,length.out)
            else x <- seq(from,to,by='days')
            i <- 1
            w <- as.POSIXlt(from)$wday
            while(from <= to){
                if ( (is.null(holidays) || !any(from==holidays)) && w != 0 && w != 6){
                    ## Day is neither a weekend day nor is it a holiday
                    ## so count as a business day
                    x[i] <- from
                    i <- i + 1
                    length.out <- length.out - 1
                }
                ## Increment/Decrement date to next day
                from <- from + 1

                ## Increment to next day, wrapping by 7
                w <- (w + 1) %% 7
            }
            x <- x[seq(len=i-1)]
        }
    }

    x
}

dateSeq.default <- dateSeq.Date
