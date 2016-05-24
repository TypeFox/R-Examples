.checkXminmax <- function (m, xminmax)
{
    new <- vapply(m, .MFclass, "")
    new <- new[names(new) %in% names(xminmax)]
    if (length(new) == 0L)
        return()
    ok <- rep(FALSE, length(new))
    names(ok) <- names(new)
    for (xvar in names(new))
        ok[xvar] <- all(m[[xvar]] >= xminmax[[xvar]][1]) && all(m[[xvar]] <= xminmax[[xvar]][2])
    if (!all(ok[!is.na(ok)])) {
        wrong <- which(!ok)
        if (sum(wrong) == 1)
            stop(gettextf("variable '%s' has values outside the permitted range",
                    names(new)[wrong]), call. = FALSE, domain = NA)
        else stop(gettextf("variables %s have values outside their permitted ranges",
                    paste(sQuote(names(new)[wrong]), collapse = ", ")),
                    call. = FALSE, domain = NA)
    }
}