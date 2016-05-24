AppendCheckDigit <- function(x, method) {
    stopifnot(is.character(x) & is.character(method) & length(method) == 1)

    FUN <- sprintf('AppendCheckDigit.%s', method)

    if (exists(FUN) && is.function(get(FUN))) {
        eval(call(FUN, x))
    } else {
        stop(sprintf('Method "%s" has not been implemented.', method))
    }
}
