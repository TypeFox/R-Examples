### ===== actuar: An R Package for Actuarial Science =====
###
### Function to discretize a continuous distribution using various
### methods.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

discretize <- function (cdf, from, to, step = 1,
                        method = c("upper", "lower", "rounding", "unbiased"),
                        lev, by = step, xlim = NULL)
{
    method <- match.arg(method)

    ## If 'cdf' is only the name of a function (say f), build a call
    ## 'f(x)'. Otherwise, check that the expression is a function call
    ## containing an 'x'. Taken from 'curve'.
    scdf <- substitute(cdf)
    if (is.name(scdf))
    {
        fcall <- paste(scdf, "(x)")
        cdf <- parse(text = fcall)
    }
    else
    {
        if (!(is.call(scdf) && match("x", all.vars(scdf), nomatch = 0)))
            stop("'cdf' must be a function or an expression containing 'x'")
        cdf <- scdf
    }

    ## If 'from' and/or 'to' are not specified, take their values in 'xlim'.
    if (missing(from))
        from <- xlim[1]
    if (missing(to))
        to <- xlim[2]

    if (method %in% c("upper", "lower"))
    {
        ## The "upper" discretization method assigns to point x =
        ## from, from + step, ..., to - step the probability mass F(x
        ## + step) - F(x).
        ##
        ## The "lower" discretization method assigns to point x = from
        ## the probability mass 0 and to x = from + step, ..., to the
        ## probability mass F(x) - F(x - step).
        ##
        ## Hence, the latter method simply has one more element than the
        ## former.
        x <- seq.int(from, to, by)
        Fx <- eval(cdf, envir = list(x = x), enclos = parent.frame())
        return(c(if(method == "lower") 0, diff(Fx)))
    }

    if (method == "rounding")
    {
        ## Rounding method assigns to point x = from the probability
        ## mass F(from + step/2) - F(from) and to point x = from +
        ## step, ..., to - step the probability mass F(x - step/2) -
        ## F(x + step/2).
        ##
        ## It is possible to make adjustments for the limits of the
        ## intervals (closed or open) for discrete distributions via
        ## 'cdf'.
        x <- c(from, seq.int(from + by/2, to - by/2, by))
        Fx <- eval(cdf, envir = list(x = x), enclos = parent.frame())
        return(diff(Fx))
    }

    if (method == "unbiased")
    {
        ## This is the matching of the first moment method. It
        ## requires a function to compute the first limited moment
        ## which should be provided in argument 'lev'. The latter is
        ## specified just like 'cdf'.
        if (missing(lev))
            stop("'lev' required with method \"unbiased\"")

        slev <- substitute(lev)
        if (is.name(slev))
        {
            fcall <- paste(slev, "(x)")
            lev <- parse(text = fcall)
        }
        else
        {
            if (!(is.call(slev) && match("x", all.vars(slev), nomatch = 0)))
                stop("'lev' must be a function or an expression containing 'x'")
            lev <- slev
        }

        ## The first limited moment must be evaluated in x = from,
        ## from + step, ..., to and the cdf in x = from and x = to
        ## only (see below).
        x <- seq.int(from, to, by)
        Ex <- eval(lev, envir = list(x = x), enclos = parent.frame())
        Fx <- eval(cdf, envir = list(x = c(from, to)), enclos = parent.frame())

        ## The probability mass in x = from is
        ##
        ##   (E[X ^ x] - E[X ^ x + step])/step + 1 - F(x).
        ##
        ## The probability mass in x = from + step, ..., to - step is
        ##
        ##   (2 * E[X ^ x] - E[X ^ x - step] - E[X ^ x + step])/step.
        ##
        ## The probability mass in x = to is
        ##
        ##   (E[X ^ x] - E[X ^ x - step])/step - 1 + F(x).
        ##
        ## See exercise 6.36 in Loss Models, 2nd edition.
        return(c(-diff(head(Ex, 2))/by + 1 - Fx[1],
                 (2 * head(Ex[-1], -1) - head(Ex, -2) - tail(Ex, -2))/by,
                 diff(tail(Ex, 2))/by - 1 + Fx[2]))
    }
}

discretise <- discretize
