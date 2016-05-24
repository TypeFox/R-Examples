### ===== actuar: An R Package for Actuarial Science =====
###
### Compute the adjustment coefficient in ruin theory, that is the
### smallest (strictly) positive root of the Lundberg equation
###
###      h(r) = E[e^(r X - r c W)] = 1,
###
### where X is the claim size random variable, W the inter-occurence
### time and c the premium rate.
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

adjCoef <- function(mgf.claim, mgf.wait = mgfexp, premium.rate, upper.bound,
                    h, reinsurance = c("none", "proportional", "excess-of-loss"),
                    from, to, n = 101)
{
    reinsurance <- match.arg(reinsurance)

    ## Sanity check
    if (missing(mgf.claim) && missing(h))
        stop("one of 'mgf.claim' or 'h' is needed")

    ## === NO REINSURANCE CASE ===
    ##
    ## Moment generating functions are unidimensional, premium rate
    ## and adjustment coefficient are both single numeric values.
    if (reinsurance == "none")
    {
        ## For each of 'mgf.claim', 'mgf.wait' and 'h' (if needed): if
        ## the expression is only the name of a function (say f),
        ## build a call 'f(x)'. Otherwise, check that the expression
        ## is a function call containing an 'x'. Taken from 'curve'
        ## and 'discretize'.
        ##
        ## NOTE: argument 'h' will be used iff 'mgf.claim' is missing,
        ## thereby giving priority to 'mgf.claim'.
        if (missing(mgf.claim))
        {
            sh <- substitute(h)
            if (is.name(sh))
            {
                fcall <- paste(sh, "(x)")
                h1 <- function(x)
                    eval(parse(text = fcall),
                         envir = list(x = x),
                         enclos = parent.frame(2))
            }
            else
            {
                if (!(is.call(sh) && match("x", all.vars(sh), nomatch = 0L)))
                    stop("'h' must be a function or an expression containing 'x'")
                h1 <- function(x)
                    eval(sh,
                         envir = list(x = x),
                         enclos = parent.frame(2))
            }
        }
        else
        {
            smgfx <- substitute(mgf.claim)
            if (is.name(smgfx))
            {
                fcall <- paste(smgfx, "(x)")
                mgfx <- parse(text = fcall)
            }
            else
            {
                if (!(is.call(smgfx) && match("x", all.vars(smgfx), nomatch = 0L)))
                    stop("'mgf.claim' must be a function or an expression containing 'x'")
                mgfx <- smgfx
            }
            smgfw <- substitute(mgf.wait)
            if (is.name(smgfw))
            {
                fcall <- paste(smgfw, "(x)")
                mgfw <- parse(text = fcall)
            }
            else
            {
                if (!(is.call(smgfw) && match("x", all.vars(smgfw), nomatch = 0L)))
                    stop("'mgf.wait' must be a function or an expression containing 'x'")
                mgfw <- smgfw
            }
            h1 <- function(x)
                eval(mgfx) * eval(mgfw, list(x = -x * premium.rate))
        }

        f1 <- function(r) (h1(r) - 1)^2

        return(optimize(f1, c(0, upper.bound - .Machine$double.eps),
                        tol = sqrt(.Machine$double.eps))$minimum)
    }

    ## === WITH REINSURANCE CASES ===
    ##
    ## Claim amount moment generating function is a function of 'x'
    ## and the retention level 'y', inter-occurence time moment
    ## generating function is a function of 'x', premium rate and
    ## adjustment coefficient are both functions of the retention
    ## level 'y'.
    ##
    ## Do same as in the no reinsurance case for each of 'mgf.claim',
    ## 'mgf.wait' and 'h' (if needed) and also 'premium'. The first
    ## must be functions of 'x' and 'y', whereas the last one is a
    ## function of 'y' only.
    if (missing(mgf.claim))
    {
        sh <- substitute(h)
        if (is.name(sh))
        {
            fcall <- paste(sh, "(x, y)")
            h2 <- function(x, y)
                eval(parse(text = fcall),
                     envir = list(x = x, y = y),
                     enclos = parent.frame(2))
        }
        else
        {
            if (!(is.call(sh) && all(match(c("x", "y"), all.vars(sh), nomatch = 0L))))
                stop("'h' must be a function or an expression containing 'x' and 'y'")
            h2 <- function(x, y)
                eval(sh,
                     envir = list(x = x, y = y),
                     enclos = parent.frame(2))
        }
    }
    else
    {
        if (!is.function(premium.rate))
            stop("'premium.rate' must be a function when using reinsurance")

        smgfx <- substitute(mgf.claim)
        if (is.name(smgfx))
        {
            fcall <- paste(smgfx, "(x, y)")
            mgfx <- parse(text = fcall)
        }
        else
        {
            if (!(is.call(smgfx) && all(match(c("x", "y"), all.vars(smgfx), nomatch = 0L))))
                stop("'mgf.claim' must be a function or an expression containing 'x' and 'y'")
            mgfx <- smgfx
        }
        smgfw <- substitute(mgf.wait)
        if (is.name(smgfw))
        {
            fcall <- paste(smgfw, "(x)")
            mgfw <- parse(text = fcall)
        }
        else
        {
            if (!(is.call(smgfw) && match("x", all.vars(smgfw), nomatch = 0L)))
                stop("'mgf.wait' must be a function or an expression containing 'x'")
            mgfw <- smgfw
        }
        spremium <- substitute(premium.rate)
        if (is.name(spremium))
        {
            fcall <- paste(spremium, "(y)")
            premium.rate <- parse(text = fcall)
        }
        else
        {
            if (!(is.call(spremium) && match("y", all.vars(spremium), nomatch = 0L)))
                stop("'premium.rate' must be a function or an expression containing 'y'")
            premium.rate <- spremium
        }

        h2 <- function(x, y)
            eval(mgfx) * eval(mgfw, list(x = -x * eval(premium.rate)))
    }

    f2 <- function(x, y) (h2(x, y) - 1)^2
    retention <- seq(from, to, length.out = n)

    ## Compute the adjustment coefficient for each retention level.
    ## The output of 'sapply' is a matrix with minima in the first
    ## line.
    ##
    ## The sapply() below passes the retention levels (argument 'y' of
    ## function 'f') to optimize(). Since the first two arguments ('f'
    ## and 'interval') of the latter function are specified, the
    ## retention levels end up in '...' and hence are considered as
    ## second argument of 'f'. *This requires R >= 2.6.0 to work since
    ## argument '...' comes much earlier in the definition of
    ## optimize().
    coef <- sapply(retention, optimize, f = f2,
                   interval = c(0, upper.bound-.Machine$double.eps),
                   tol = sqrt(.Machine$double.eps))[1L, ]

    ## Make a function from the (retention, coefficient) pairs
    ## computed above, joining the points by straight line segments.
    FUN <- approxfun(retention, coef, rule = 2, method = "linear")

    comment(FUN) <- paste(toupper(substring(reinsurance, 1L, 1L)),
                          substring(reinsurance, 2L),
                          " reinsurance",
                          sep = "", collapse = "")
    class(FUN) <- c("adjCoef", class(FUN))
    attr(FUN, "call") <- sys.call()
    FUN
}

plot.adjCoef <- function(x, xlab = "x", ylab = "R(x)",
                         main = "Adjustment Coefficient",
                         sub = comment(x), type = "l", add = FALSE, ...)
{
    xx <- eval(expression(x), envir = environment(x))
    yy <- eval(expression(y), envir = environment(x))
    if (add)
        lines(xx, yy, ..., main = main, xlab = xlab, ylab = ylab,
              type = type)
    else
        plot(xx, yy, ..., main = main, xlab = xlab, ylab = ylab,
             type = type)
    mtext(sub, line = 0.5)
}
