isSuperSet <- function(x, y){
    ## x and y is a vector representing a primeImplicant, e.g. x <- c(1,-9,-9); y <- c(1,1,1)
    ## return TRUE if x is a superset of y
    idx <- !is.dontcare(x)
    equal <- x[idx]==y[idx]
    dontcare <- is.dontcare(y[idx])
    if (all(equal | dontcare)) {
        TRUE
    } else {
        FALSE
    }
}
