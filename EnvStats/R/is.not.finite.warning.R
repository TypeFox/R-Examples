is.not.finite.warning <-
function (x, name = deparse(substitute(x))) 
{
    if (!is.numeric(x)) {
        warning(paste(name, "is not a numeric dataset"))
    }
    else {
        n.na <- sum(is.na(x))
        n.nan <- sum(is.nan(x))
        n.inf <- sum(is.infinite(x))
        n.bad <- n.na + n.inf
        if (n.bad == 0) 
            return()
        n <- c(n.na - n.nan, n.nan, n.inf)
        types <- c("NA's", "NaN's", "+-Inf's")
        msg0 <- paste("There were", n.bad, "nonfinite values in", 
            name, ":")
        msg1 <- paste(n[n > 0], types[n > 0], collapse = ", ")
        warning(paste(msg0, msg1))
    }
}
