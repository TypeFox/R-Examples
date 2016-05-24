##
## Below are versions of the 'c' function for POSIXt derived classes that
## do not purposefully throw away the timezone, like R does. Our
## versions assume that all the passed in POSIXt objects are in the
## same timezone.
##
## They avoid the following R gotchas:
##
## ## Works fine in Date only
## > do.call("c", lapply(0:4, function(i) dateParse("2008/12/25") + i))
## [1] "2008-12-25" "2008-12-26" "2008-12-27" "2008-12-28" "2008-12-29"
##
## ## 'as.POSIXlt' goes to UST. 'c' goes to local time zone??????????
## > do.call("c", lapply(0:4, function(i) as.POSIXlt(dateParse("2008/12/25") + i)))
## [1] "2008-12-24 17:00:00 MST" "2008-12-25 17:00:00 MST"
## [3] "2008-12-26 17:00:00 MST" "2008-12-27 17:00:00 MST"
## [5] "2008-12-28 17:00:00 MST"
##
## ## Use our 'combine' instead. They preserve the timezone.
## > do.call("combine", lapply(0:4, function(i) as.POSIXlt(dateParse("2008/12/25") + i)))
## [1] "2008-12-25 UTC" "2008-12-26 UTC" "2008-12-27 UTC" "2008-12-28 UTC"
## [5] "2008-12-29 UTC"
##

combine <- function (object, ...) UseMethod("combine")

combine.POSIXct <- function(..., recursive = FALSE)
{
    l <- list(...)
    saveAttr <- attributes(l[[1]])
    res <- c(unlist(lapply(list(...), unclass)))
    attributes(res) <- saveAttr
    res
}

combine.POSIXlt <- function(..., recursive = FALSE)
{
    as.POSIXlt(do.call("combine", lapply(list(...), as.POSIXct)))
}
