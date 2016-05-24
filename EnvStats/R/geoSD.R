geoSD <-
function (x, na.rm = FALSE, sqrt.unbiased = TRUE) 
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
    else {
        sd.log <- sd(log(x))
        if (!sqrt.unbiased) {
            n <- length(x)
            sd.log <- sqrt((n - 1)/n) * sd.log
        }
    }
    exp(sd.log)
}
