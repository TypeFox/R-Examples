########################
##  tclust.opt.restr  ##  ##!!  not to be released yet - still beta  !!##
########################

.absfact <- function (x, y)
{
    x <- abs (x)
    y <- abs (y)
    if (x > y)
        x / y 
    else
        y / x
}

.tclust.opt.restr <- function (..., k, alpha, max.restr.fact = 1e3, trace = 0)
{    
    max.rf <- max.restr.fact
    init.rf <- 1
    mult.rf <- 1.1
    inc.rf <- 1.5
    max.inc <- 4

    pass.trace <- ifelse (trace > 2, trace - 2, 0)
    cur.rf <- init.rf

    repeat
    {
        if (trace >= 2)
            cat ("checking restriction factor", round (cur.rf, 2), "\n")
        r <- tclust (..., k = k, alpha = alpha, restr.fact = cur.rf, trace = max (trace - 2, 0), warnings = FALSE)

        if (r$unrestr.fact < cur.rf)
            break
        if (cur.rf > max.rf)
            stop ("could not find an appropriate restriction factor. Try increasing max.rf")
        cur.rf <- min (cur.rf * max.inc, r$unrestr.fact) * inc.rf
    }

    lrf <- round (r$unrestr.fact * mult.rf, 2)

    cur.rf <- max.rf

    max.it <- 5

    last.rf <- NULL
    for (i in 1:max.it)    ##    after ~5 iterations we should have found the factor
    {
        if (trace >= 2)    
            cat ("checking restriction factor", round (cur.rf, 2), "\n")

        r <- tclust (..., k = k, alpha = alpha, restr.fact = cur.rf, trace = max (trace - 2, 0), warnings = FALSE)

        if (.absfact (cur.rf, r$unrestr.fact) <= inc.rf)
            break

        if (!is.null (last.rf) &&
            .absfact (last.rf, r$unrestr.fact) <= inc.rf)
            break

        last.rf <- r$unrestr.fact
        cur.rf <- min (cur.rf * max.inc, r$unrestr.fact) * inc.rf
    }

    if (cur.rf > max.rf || i == max.it)
        .stop.unstable (alpha, k)

    urf <- round (r$unrestr.fact * mult.rf, 2)

    r$conv <- .absfact (urf, lrf) <= inc.rf

    if (!r$conv)
    {
        warning (paste ("Search for restriction factor did not converge (",
                 round(lrf, 2), " vs. ", round (urf, 2), 
                 ")\n  Use \"DiscrFact\" to ensure the result's integrity.", 
                 sep = "")
                 )
    }

    r$restr.fact <- ceiling (urf)

    if (trace >= 1)    
        cat ("found restriction factor", r$restr.fact, "\n")

    return (r)
}

.stop.unstable <- function (alpha, k)
{
      stop (paste ("The chosen parameters alpha =", round (alpha, 2), 
            "and k =", k, "appear to yield spurious clusters.",
            "\n  Try reducing one of these parameters."
#            ," in order to get a proper selection of \"restr.fact\""
            ))
}
