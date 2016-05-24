
returnMessage <- function(x, ...)
    UseMethod("returnMessage")

returnMessage.default <- function(x, ...)
    x$returnMessage

returnMessage.maxim <- function(x, ...)
    x$message

returnMessage.maxLik <- function(x, ...)
    x$message
