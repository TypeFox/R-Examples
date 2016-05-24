famfun <- function(fam, deriv, theta) {
    setfam(list(fam))
    result <- .C("aster_family",
        index = as.integer(1),
        deriv = as.integer(deriv),
        theta = as.double(theta),
        value = double(1),
        PACKAGE = "aster")$value
    clearfam()
    return(result)
}
