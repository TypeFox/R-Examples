### ===== actuar: An R Package for Actuarial Science =====
###
### Create modified density and modified cumulative distribution
### function for data with deductible, limit, coinsurance and
### inflation.
###
### See Chapter 5 of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

coverage <- function(pdf, cdf, deductible = 0, franchise = FALSE,
                     limit = Inf, coinsurance = 1, inflation = 0,
                     per.loss = FALSE)
{
    Call <- match.call()

    ## First determine if the cdf is needed or not. It is needed when
    ## there is a deductible or a limit and, of course, if the output
    ## function should compute the cdf.
    is.cdf <- missing(pdf) || is.null(pdf) # return cdf?
    has.limit <- limit < Inf               # used often
    needs.cdf <- any(deductible > 0, has.limit, is.cdf) # cdf needed?

    ## Sanity check of arguments
    if (any(deductible < 0, limit < 0, coinsurance < 0, inflation < 0))
        stop("coverage modifications must be positive")
    if (limit <= deductible)
      stop("deductible must be smaller than the limit")
    if (coinsurance > 1)
        stop("coinsurance must be between 0 and 1")
    if (missing(cdf) & needs.cdf)
        stop("'cdf' must be supplied")

    ## Quantites often used. Leave as expressions for the output
    ## function to preserve accuracy.
    r <- 1 + inflation
    d <- if (inflation) substitute(d/r, list(d = deductible, r = r)) else deductible
    u <- if (inflation) substitute(u/r, list(u = limit, r = r)) else limit

    ## The modified cdf or pdf are usually defined in branches. To
    ## avoid using nested ifelse(), we will rather rely on sets of
    ## expressions to make the required calculations for each branch
    ## separately. This is actually much faster.
    ##
    ## The output function will have varying number of said
    ## expressions depending on the case that is dealt with. We will
    ## build the body of the output function piece by piece as we go
    ## along.
    e <- expression(Call <- match.call())

    ## One main discriminating factor is whether the cdf is needed for
    ## the output function of not.
    if (needs.cdf)
    {
        ## Get argument list of 'cdf' to transfert them to the output
        ## function.
        argv <- formals(cdf)            # arguments as list
        argn <- names(argv)             # arguments names as strings

        ## Remember if argument 'lower.tail' is available, so we can
        ## use it later. Then, drop unsupported arguments 'lower.tail'
        ## and 'log.p'.
        has.lower <- "lower.tail" %in% argn
        argn <- setdiff(argn, c("lower.tail", "log.p"))

        ## Calculations in the output function are done by evaluating
        ## function calls built upon the call to the output function
        ## itself. This is convenient as we do not need to fiddle with
        ## the values of the formal arguments.
        if (is.cdf)                     # output function computes F(y)
        {
            ## The output function will have the same formal arguments
            ## as the cdf. Object 'x' holds the symbol of the first
            ## argument.
            argsFUN <- argv[argn]       # arguments of output function
            x <- as.name(argn[1])       # symbol

            ## Calls needed in this case:
            ## 1. one to compute F(y);
            ## 2. one to compute F(d) if there is a deductible;
            ## 3. one to compute 1 - F(d) for the per payment cases.
            ## Never need to compute F(u).
            ##
            ## If 'lower.tail' is available in 'cdf', then set it to
            ## FALSE to compute S(d) = 1 - F(d) more accurately.
            e <- c(e,
                   quote(F <- Call),        # 1.
                   substitute(F[[1L]] <- as.name(fun),
                              list(fun = as.character(Call$cdf))))
            if (deductible)
            {
                e <- c(e,
                       quote(Fd <- F),      # 2.
                       substitute(Fd[[2L]] <- a, list(a = d)))
                if (!per.loss & has.lower)
                    e <- c(e,
                           quote(Sd <- Fd), # 3.
                           quote(Sd$lower.tail <- FALSE))
            }
        }
        else                            # output function computes f(y)
        {
            ## When there is a limit, we will need to compute 1 - F(u)
            ## as is or using 'lower.tail = FALSE' to improve
            ## accuracy. For clarity in the output function, we will
            ## use Fu as object name in the first case and Su in the
            ## second.
            if (has.limit)
            {
                if (has.lower)
                {
                    Fu.name <- as.name("Su")
                    Su.quote <- quote(eval.parent(Su))
                }
                else
                {
                    Fu.name <- as.name("Fu")
                    Su.quote <- quote((1 - eval.parent(Fu)))
                }
            }

            ## Calls needed in this case:
            ## 1. one to compute F(d) if there is a deductible for the
            ##    per loss cases;
            ## 2. one to compute 1 - F(d) if there is a deductible for the
            ##    per payment cases;
            ## 3. one to compute 1 - F(u) when there is a limit.
            ## No function to compute F(y) needed.
            ##
            ## If 'lower.tail' is available in 'cdf', then set it to
            ## FALSE to compute S(d) = 1 - F(d) and S(u) = 1 - F(u)
            ## more accurately.
            if (deductible)             # f(y) with deductible
            {
                Fd.name <- as.name(if (!per.loss & has.lower) "Sd" else "Fd")
                e <- c(e,
                       substitute(G <- Call, list(G = Fd.name)),         # 1. or 2.
                       if (!per.loss & has.lower) quote(Sd$lower.tail <- FALSE),
                       substitute(G[[1L]] <- as.name(fun),
                                  list(G = Fd.name, fun = as.character(Call$cdf))),
                       substitute(names(G)[2L] <- q, list(G = Fd.name, q = argn[1])),
                       substitute(G[[2L]] <- a, list(G = Fd.name, a = d)))
                if (has.limit)
                    e <- c(e,
                           substitute(H <- G, list(H = Fu.name, G = Fd.name)), # 3.
                           if (per.loss & has.lower) quote(Su$lower.tail <- FALSE),
                           substitute(H[[2L]] <- a, list(H = Fu.name, a = u)))
            }
            else                        # f(y) with limit only
            {
                ## Since 'needs.cdf == TRUE', then this case
                ## necessarily has 'limit < Inf'. Only call needed is
                ## one to compute 1 - F(u).
                e <- c(e,
                       substitute(G <- Call, list(G = Fu.name)),
                       if (has.lower) quote(Su$lower.tail <- FALSE),
                       substitute(G[[1L]] <- as.name(fun),
                                  list(G = Fu.name, fun = as.character(Call$cdf))),
                       substitute(names(G)[2L] <- q, list(G = Fu.name, q = argn[1])),
                       substitute(G[[2L]] <- a, list(G = Fu.name, a = u)))
            }
        }
    }

    ## Repeat same steps as above for case needing the pdf. The output
    ## function is a pdf and in this case the arguments of the output
    ## function are those of 'pdf'.
    if (!is.cdf)
    {
        argv <- formals(pdf)                # arguments as list
        argn <- setdiff(names(argv), "log") # drop argument 'log'
        argsFUN <- argv[argn]           # arguments of output function
        x <- as.name(argn[1])           # symbol
        e <- c(e,
               quote(f <- Call),
               substitute(f[[1L]] <- as.name(fun),
                          list(fun = as.character(Call$pdf))))
    }

    ## Build the value at which the underlying pdf/cdf will be called
    ## for non special case values of 'x'. We need to index 'x' to
    ## only compute for the correct values of a given branch.
    x.mod <- as.call(c(as.name("["), x, as.name("w")))
    if (coinsurance < 1)
        x.mod <- substitute(x/alpha, list(x = x.mod, alpha = coinsurance))
    if (deductible & !franchise)
        x.mod <- substitute(x + d, list(x = x.mod, d = deductible))
    if (inflation)
        x.mod <- substitute((x)/r, list(x = x.mod, r = r))

    ## Each pdf/cdf is defined in three branches. Define the
    ## boundaries and conditions for the first two branches. Those for
    ## the third branch are defined further down.
    if (franchise)
    {
        bound1 <- coinsurance * deductible
        bound2 <- coinsurance * limit
        cond1 <- if (is.cdf)
            substitute(0 <= x & x <= b1, list(x = x, b1 = bound1))
        else
            quote(x == 0)
        cond2 <- substitute(b1 < x & x < b2,
                            list(x = x, b1 = bound1, b2 = bound2))
    }
    else
    {
        bound1 <- 0
        bound2 <- coinsurance * (limit - deductible)
        cond1 <- substitute(x == 0, list(x = x))
        cond2 <- substitute(0 < x & x < b, list(x = x, b = bound2))
    }

    ## Initialization of the results vector in the output function
    ## with 0s.
    e <- c(e,
           substitute(res <- numeric(length(x)), list(x = x)))

    ## Definition of the output function for the first branch. There
    ## is a computation to make only if there is a deductible with the
    ## payment per loss random variable. For all other cases, the
    ## value in the first branch is 0 and we rely on the
    ## initialization with numeric() done at the previous step.
    if (per.loss & deductible)
        e <- c(e,
               substitute(res[which(cond1)] <- eval.parent(Fd),
                          list(cond1 = cond1)))

    ## Definition of the output function for the second and third
    ## branches. The 'is.cdf = TRUE' and 'is.cdf = FALSE' cases must
    ## be treated separately.
    if (is.cdf)
    {
        cond3 <- substitute(x >= b, list(x = x, b = bound2))
        f2 <- quote(eval.parent(F))
        if (!per.loss & deductible)
            f2 <- if (has.lower)
                substitute((f - F)/S,
                           list(f = f2,
                                F = quote(eval.parent(Fd)),
                                S = quote(eval.parent(Sd))))
            else
                substitute((f - F)/S,
                           list(f = f2,
                                F = quote((p <- eval.parent(Fd))),
                                S = quote((1 - p))))
        e <- c(e,
               substitute(w <- which(cond), list(cond = cond2)),
               substitute(F[[2L]] <- x, list(x = x.mod)),
               substitute(res[w] <- f, list(f = f2)),
               if (has.limit) substitute(res[cond] <- 1, list(cond = cond3)))
    }
    else
    {
        cond3 <- substitute(x == b, list(x = x, b = bound2))
        f2 <- quote(eval.parent(f))
        if (has.limit) f3 <- Su.quote
        if (!per.loss & deductible)
        {
            if (has.limit)
            {
                f2 <- if (has.lower)
                    substitute(f/(p <- S),
                               list(f = f2, S = quote(eval.parent(Sd))))
                else
                    substitute(f/(p <- S),
                               list(f = f2, S = quote(1 - eval.parent(Fd))))
                f3 <- substitute(f/p, list(f = f3))
            }
            else
                f2 <- if (has.lower)
                    substitute(f/S,
                               list(f = f2, S = quote(eval.parent(Sd))))
                else
                    substitute(f/S,
                               list(f = f2, S = quote((1 - eval.parent(Fd)))))
        }
        if (inflation | coinsurance < 1)
            f2 <- substitute(f/k, list(f = f2, k = coinsurance * r))
        e <- c(e,
               substitute(w <- which(cond), list(cond = cond2)),
               substitute(f[[2L]] <- x, list(x = x.mod)),
               substitute(res[w] <- f, list(f = f2)),
               if (has.limit) substitute(res[cond] <- f, list(cond = cond3, f = f3)))
    }

    ## Last expression of the output function.
    e <- c(e, quote(res))

    ## Wrap up the output function.
    FUN <- function() {}
    body(FUN) <- as.call(c(as.name("{"), e)) # taken from help(body)
    formals(FUN) <- argsFUN             # set arguments
    environment(FUN) <- new.env()       # new, empty environment
    FUN
}
