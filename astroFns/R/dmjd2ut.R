dmjd2ut <-
    function(dmjd, tz='UTC') {

# Decimal Modified Julian Date to UT calculation
# tz is an optional time zone string (e.g. in US EST5EDT or EST, etc.)
# number of digits in seconds display, n, is controlled by options('digits.secs'=n)
# A. Harris, U. Maryland Astronomy, 3/17/2008 to 9/19/2012

    # Convert to Julian date and compute year, month, and day
    jd <- trunc(dmjd) + 2400000.5  # here jd is 0h
    ymd <- jd2ymd(jd)
    ymd <- as.numeric(unlist(strsplit(as.character(ymd), '[-: ]')))
    len <- length(ymd)
    yr <- ymd[seq(1, len, by=3)]
    mo <- ymd[seq(2, len, by=3)]
    dy <- ymd[seq(3, len, by=3)]
    # Work out hours, minutes and seconds; round seconds for good display
    dayfrac <- dmjd%%1
    hr <- dayfrac*24
    min <- (hr%%1)*60
    ds <- getOption('digits.secs')
    if(is.null(ds)) ds <- 0
    sec <- round((min%%1)*60, ds)

    # Work out date in UTC
    out <- ISOdatetime(yr, mo, dy, trunc(hr), trunc(min), sec, 'UTC')
    # Then change to appropriate time zone
    attr(out, 'tzone') <- tz
    out
}

