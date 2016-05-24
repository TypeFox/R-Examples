`confint.cusp` <-
function (object, parm, level = 0.95, ...) 
{
    cf <- coef(object)
    pnames <- names(cf)
    v <- vcov(object)
    pnames <- colnames(v)
    if (missing(parm)) 
        parm <- seq_along(pnames)
    else if (is.character(parm)) 
        parm <- match(parm, pnames, nomatch = 0)
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- format.perc(a, 3)
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(parm), 2), dimnames = list(pnames[parm], 
        pct))
    ses <- sqrt(diag(v))[parm]
    ci[] <- cf[pnames[parm]] + ses %o% fac
    ci
}

