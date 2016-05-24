.getXminmax <- function (Terms, m)
{
    deparse2 <- function(x) paste(deparse(x, width.cutoff = 500L),
        collapse = " ")
    xvars <- sapply(attr(Terms, "variables"), deparse2)[-1L]
    if ((yvar <- attr(Terms, "response")) > 0)
        xvars <- xvars[-yvar]
    if (length(xvars)) {
        xminmax <- lapply(m[xvars], function(x) if (is.numeric(x))
            c(min(x), max(x))
        else NULL)
        xminmax[!vapply(xminmax, is.null, NA)]
    }
    else NULL
}