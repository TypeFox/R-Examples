## Returns return (error) code

returnCode <- function(x, ...)
    UseMethod("returnCode")

returnCode.default <- function(x, ...)
    x$returnCode

returnCode.maxLik <- function(x, ...)
    x$code
