##
## Unfortunately, 'format.POSIXct' is not compatible between
## Linux and Windows R. Here are a few for the difference I have found:
## 1) Windows version does not recognize "%y" format.
## 2) Windows version does not recognize width arguments like
##    "%02d".
## 3) Windows and Linux does not agree on the meaning of "%Y".
##    Under Windows it means "%04Y"; under Linux it prints with
##    minimal width.
##
## This function formats any kind of data objects to character
## strings with the default format "%02m-%02d-%04Y" under
## both Linux and Windows.
##

dateFormat.env <- new.env()
dateFormat.env$initialized <- FALSE

dateFormat.config <- function() {
    # figure out the capabilities of the date formatting on this system
    # Use the year 0045 as a test case.  Some systems that accept a %04Y
    # format will format %Y 1991 as '1991', but %Y 0045 as '45'.
    if (strftime(as.POSIXct('0045-01-02'), '%Y-%m-%d') == '0045-01-02')
        dateFormat.env$default.format <- '%Y-%m-%d'
    else if (strftime(as.POSIXct('0045-01-02'), '%04Y-%m-%d') == '0045-01-02')
        dateFormat.env$default.format <- '%04Y-%m-%d'
    else if (strftime(as.POSIXct('0045-01-02'), '%04Y-%02m-%02d') == '0045-01-02')
        dateFormat.env$default.format <- '%04Y-%02m-%02d'
    else if (strftime(as.POSIXct('0045-01-02'), '%Y-%02m-%02d') == '0045-01-02')
        dateFormat.env$default.format <- '%Y-%02m-%02d'
    else
        dateFormat.env$default.format <- '%Y-%m-%d'

    if (strftime(as.POSIXct('0045-01-02'), '%04Y')=='0045')
        dateFormat.env$f04Y.ok <- TRUE
    else
        dateFormat.env$f04Y.ok <- FALSE

    if (strftime(as.POSIXct('1991-01-02'), '%02m-%02d')=='01-02')
        dateFormat.env$f02md.ok <- TRUE
    else
        dateFormat.env$f02md.ok <- FALSE
}

dateFormat <- function(date, format = NULL)
{
    if (!dateFormat.env$initialized)
        dateFormat.config()
    if (is.null(format))
        format <- dateFormat.env$default.format
    if (!dateFormat.env$f04Y.ok && regexpr('%04Y', format, fixed=TRUE)>0)
        format <- gsub('%04Y', '%Y', format, fixed=TRUE)
    if (!dateFormat.env$f02md.ok && regexpr('%02', format, fixed=TRUE)>0)
        format <- gsub('%02', '%', format, fixed=TRUE)

    if (is.character(date))
        date <- dateParse(date)

    if (inherits(date, "dates"))
        date <- as.POSIXct(date)
    if (!(inherits(date, "Date") || is(date, "POSIXt")))
        stop("unknown date class: '", class(date), "'")

    # Intervene and insert formats %Q and %C, but correctly ignore %%
    if (any(i <- regexpr('%+Q', format))>=1 && any(attr(i, 'match.length')%%2==0)) {
        if (length(format) < length(date)) {
            format <- rep(format, len=length(date))
            i <- structure(rep(i, len=length(date)), match.length=rep(attr(i, "match.length"), len=length(date)))
        }
        j <- (attr(i, 'match.length')%%2)==0
        if (any(attr(i, 'match.length')[j]>2))
            i[j] <- i[j] + attr(i, 'match.length')[j]-2
        substring(format[j], i[j], i[j]+2) <- quarters(date)
    }
    if (all(regexpr('%', format)<=0))
        return(rep(format, len=length(date)))
    if (.Platform$OS.type=='windows' && any((i <- regexpr('%+y', format))>=1) && any(attr(i, 'match.length')%%2==0)
        && as.numeric(format(min(date), '%Y')) < 1900) {
        # '%y' can be buggy on windows for dates < 1900
        # e.g., format(as.Date('1899-01-01'), '%y') returns "0/" on a Win XP 64 bit OS in 2012
        if (length(format) < length(date)) {
            format <- rep(format, len=length(date))
            i <- structure(rep(i, len=length(date)), match.length=rep(attr(i, "match.length"), len=length(date)))
        }
        j <- (attr(i, 'match.length')%%2)==0
        if (any(attr(i, 'match.length')[j]>2))
            i[j] <- i[j] + attr(i, 'match.length')[j]-2
        substring(format[j], i[j], i[j]+2) <- sprintf('%02g', as.numeric(format(date, '%Y')) %% 100)
    }
    if (all(regexpr('%', format)<=0))
        return(rep(format, len=length(date)))
    if (any((i <- regexpr('%+C', format))>=1) && any(attr(i, 'match.length')%%2==0)) {
        if (length(format) < length(date)) {
            format <- rep(format, len=length(date))
            i <- structure(rep(i, len=length(date)), match.length=rep(attr(i, "match.length"), len=length(date)))
        }
        j <- (attr(i, 'match.length')%%2)==0
        if (any(attr(i, 'match.length')[j]>2))
            i[j] <- i[j] + attr(i, 'match.length')[j]-2
        substring(format[j], i[j], i[j]+2) <- sprintf('%02g', as.numeric(format(date, '%Y')) %/% 100)
    }
    if (all(regexpr('%', format)<=0))
        return(rep(format, len=length(date)))
    if (length(format) > 1)
        return(mapply(strftime, date, format))
    else
        return(strftime(date, format))
}
