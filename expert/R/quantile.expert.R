### ===== experts =====
###
### Quantiles for objects of class "expert"
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
###          Mathieu Pigeon <mathieu.pigeon.3@ulaval.ca>

quantile.expert <- function(x, probs = seq(0, 1, 0.25),
                            smooth = FALSE, names = TRUE, ...)
{
    y <- c(0, cumsum(x$probs))
    x <- x$breaks
    ind <- sapply(probs, function(q) match(TRUE, y >= q))

    ## Create the inverse function of either the cdf or the ogive.
    fun <-
        if (smooth)                     # ogive
            approxfun(y, x, yleft = min(x), yright = max(x),
                      method = "linear", ties = "ordered")
        else                            # cdf
            approxfun(y, x, yleft = min(x), yright = max(x),
                      method = "constant", ties = "ordered")

    ## Quantiles
    res <- fun(probs)

    if (names)
    {
        dig <- max(2, getOption("digits"))
        names(res) <- formatC(paste(100 * probs, "%", sep = ""),
                              format = "fg", wid = 1, digits = dig)
    }
    res
}
