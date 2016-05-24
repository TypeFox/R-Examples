parseDollars <- function(x, pattern='\\$|,', replacement='', ...){
    X <- gsub(pattern, replacement, x, ...)
    as.numeric(X)
}

