iqr <-
function (x, na.rm = FALSE) 
{
    wna <- which.na(x)
    if (length(wna)) {
        if (na.rm) 
            x <- x[-wna]
        else return(NA)
    }
    ret.val <- diff(quantile(x, probs = c(0.25, 0.75)))
    names(ret.val) <- NULL
    ret.val
}
