###
### $Id: strcmp.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Compare strings.
###


##-----------------------------------------------------------------------------
strcmp <- function(S, T) {
    if (!is.character(S)) {
        stop(sprintf("argument %s must be character", sQuote("S")))
    }
    if (!is.character(T)) {
        stop(sprintf("argument %s must be character", sQuote("T")))
    }

    return(if (length(S) == length(T)) {
               all(S == T)
           } else {
               FALSE
           })
}

