## Unexported Utility Functions


#' Finds the first value of a vector that is less than a value.
#' 
#' None
#' 
#' @param x Vector to serach
#' @param val Threshhold
#' @return Gives the index of the first value of x that is <= val.
#' returns -1 if no value satisfies the condition
indexFinder <- function(x, val) {                 
    idx <- (1:length(x))[x<= val][1]
    if (is.na(idx)) {idx <- -1}
    idx
}
