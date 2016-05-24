years <- function(x) as.POSIXlt(x)$year + 1900

emptyDate <- function() .dateParse.origin[0]

is.wholenumber <- function(x) any(floor(x) == x)

