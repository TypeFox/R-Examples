parseCommas <- function(x, pattern='\\$|,', replacement='',
                                acceptableErrorRate=0, ... ){
    UseMethod('parseCommas')
}

parseCommas.default <- function(x, pattern='\\$|,', replacement='',
                                acceptableErrorRate=0, ... ){
    if(is.factor(x) | is.character(x)){
        xc <- as.character(x)
        x. <- parseDollars(xc, pattern, replacement, ...)
        good <- mean(is.na(x.[!is.na(x)]))
        if(good <= acceptableErrorRate) {
            return(x.)
        } else return(x)
    }
    x
}

parseCommas.data.frame <- function(x, pattern='\\$|,',
              replacement='', acceptableErrorRate=0, ... ){
    X <- lapply(x, parseCommas)
    as.data.frame(X)
}
