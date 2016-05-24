geoMean <-
function (x, na.rm = FALSE) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    wna <- which.na(x)
    if (length(wna)) {
        if (na.rm) 
            x <- x[-wna]
        else return(NA)
    }
    if (any(x <= 0)) {
        warning("Non-positive values in 'x'")
        return(NA)
    }
    else return(exp(mean(log(x))))
}
