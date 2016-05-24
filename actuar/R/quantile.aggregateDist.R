### ===== actuar: An R Package for Actuarial Science =====
###
### Quantiles for objects of class 'aggregateDist'
###
### AUTHORS:  Louis-Philippe Pouliot,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

quantile.aggregateDist <-
    function(x, probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.995),
             smooth = FALSE, names = TRUE, ...)
{
    chkDots(...)                        # method does not use '...'
    label <- comment(x)

    ## The Normal and Normal Power approximations are the only
    ## continuous distributions of class 'aggregateDist'. They are
    ## therefore treated differently, using the 'base' quantile
    ## function qnorm().
    if (label == "Normal approximation")
        res <- qnorm(probs, get("mean", environment(x)),
                     sqrt(get("variance", environment(x))))
    else if (label == "Normal Power approximation")
    {
        m <- get("mean", envir = environment(x))
        sd <- sqrt(get("variance", envir = environment(x)))
        sk <- get("skewness", envir = environment(x))

        ## Calling qnorm() and inverting the Normal Power 'standardization'
        q <- qnorm(probs)
        res <- ifelse(probs <= 0.5, NA, m + sd * (q + sk * (q^2 - 1)/6))
    }
    else
    {
        ## An empirical and discrete approach is used for
        ## 'aggregateDist' objects obtained from methods other than
        ## Normal and Normal Power.
        y <- get("y", environment(x))
        x <- get("x", environment(x))

        ## Create the inverse function of either the cdf or the ogive.
        fun <-
            if (smooth)                     # ogive
                approxfun(y, x, yleft = 0, yright = max(x),
                          method = "linear", ties = "ordered")
            else                            # cdf
                approxfun(y, x, yleft = 0, yright = max(x),
                          method = "constant", f = 1, ties = "ordered")

        ## Quantiles
        res <- fun(probs)
    }

    if (names)
    {
        dig <- max(2, getOption("digits"))
        names(res) <- formatC(paste(100 * probs, "%", sep = ""),
                              format = "fg", width = 1, digits = dig)
    }
    res
}
