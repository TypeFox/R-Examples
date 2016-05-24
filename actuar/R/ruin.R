### ===== actuar: An R Package for Actuarial Science =====
###
### Calulation of infinite time ruin probabilities in the models of
### Cramer-Lundberg and Sparre Andersen. In general, one has
###
###   psi(u) = pphtype(u, pi, Q, lower.tail = FALSE)
###
### for definitions of pi and Q that depend on the severity and waiting
### times models. An explicit solution exists when both severity and
### waiting times are exponential (no mixture).
###
### _Combinations_ of exponentials as defined in Dufresne & Gerber
### (1988) are NOT supported.
###
### References:
###
### Dufresne, F. and Gerber, H. U. (1988), "Three methods to calculate
### the probability of ruin", Astin Bulletin 19, p. 71-90.
###
### Asmussen, S. and Rolski, T. (1991), "Computational methods in risk
### theory: A matrix-algorithmic approach", Insurance: Mathematics and
### Economics 10, p. 259-274.
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

ruin <- function(claims = c("exponential", "Erlang", "phase-type"), par.claims,
                 wait = c("exponential", "Erlang", "phase-type"), par.wait,
                 premium.rate = 1, tol = sqrt(.Machine$double.eps), maxit = 200,
                 echo = FALSE)
{
    ## Sanity checks
    if (missing(par.claims) || !is.list(par.claims))
        stop("'par.claims' must be a named list")
    if (missing(par.wait) || !is.list(par.wait))
        stop("'par.wait' must be a named list")

    claims <- match.arg(claims)
    wait <- match.arg(wait)

    ## ==============================================
    ##  Extraction of the claims severity parameters
    choices <-
        switch(claims,
               exponential = c("rate", "weights"),
               Erlang = c("shape", "rate", "scale", "weights"),
               "phase-type" = c("prob", "rates"))
    i <- pmatch(names(par.claims), choices, nomatch = 0L, duplicates.ok = TRUE)
    if (all(i == 0L))
        stop(gettextf("parameters %s missing in 'par.claims'",
                      paste(dQuote(choices), collapse = ", ")))
    par.claims <- par.claims[i > 0L]    # keep relevant components
    p <- choices[i[i > 0L]]             # keep relevant names
    names(par.claims) <- p              # use full names

    if (claims == "exponential")
    {
        if ("rate" %in% p)
            rate <- par.claims$rate
        else
            stop("parameter \"rate\" missing in 'par.claims'")
        n <- length(rate)

        if ("weights" %in% p && n > 1L)
            prob <- rep(par.claims$weights, length.out = n)
        else if (n == 1L)
            prob <- 1
        else
            stop("parameter \"weights\" missing in 'par.claims'")

        rates <- diag(-rate, n)
    }
    else if (claims == "Erlang")
    {
        if ("shape" %in% p)
            shape <- par.claims$shape
        else
            stop("parameter \"shape\" missing in 'par.claims'")
        if ("rate" %in% p)
            rate <- par.claims$rate
        else if ("scale" %in% p)
            rate <- 1 / par.claims$scale
        else
            stop("parameter \"rate\" or \"scale\" missing in 'par.claims'")
        if (length(shape) < length(rate))
            shape <- rep(shape, length.out = length(rate))
        else
            rate <-  rep(rate, length.out = length(shape))
        n <- sum(shape)

        if ("weights" %in% p && length(shape) > 1L)
        {
            prob <- numeric(n)
            prob[cumsum(c(1, head(shape, -1)))] <- par.claims$weights
        }
        else if (length(shape) == 1L)
            prob <- c(1, rep(0, n - 1))
        else
            stop("parameter \"weights\" missing in 'par.claims'")

        rates <- diag(rep(-rate, shape), n)
        if (n > 1 && shape > 1)
        {
            tmp <- -head(diag(rates), -1L)
            tmp[cumsum(head(shape, -1L))] <- 0 # insert 0s in "ll corners"
            rates[cbind(seq_len(n - 1), seq(2, len = n - 1))] <- tmp
        }
    }
    else                                # claims == "phase-type"
    {
        if ("prob" %in% p)
            prob <- par.claims$prob
        else
            stop("parameter \"prob\" missing in 'par.claims'")
        if ("rates" %in% p)
            rates <- par.claims$rates
        else
            stop("parameter \"rates\" missing in 'par.claims'")
        n <- length(prob)
        if (!(is.matrix(rates) && nrow(rates) == n))
            stop("invalid parameters in 'par.claims'")
    }
    ## ==============================================

    ## =================================================
    ##  Extraction of the interarrival times parameters
    choices <-
        switch(wait,
               exponential = c("rate", "weights"),
               Erlang = c("shape", "rate", "scale", "weights"),
               "phase-type" = c("prob", "rates"))
    i <- pmatch(names(par.wait), choices, nomatch = 0L, duplicates.ok = TRUE)
    if (all(i == 0L))
        stop(gettextf("parameters %s missing in 'par.wait'",
                      paste(dQuote(choices), collapse = ", ")))
    par.wait <- par.wait[i > 0L]        # keep relevant components
    p <- choices[i[i > 0L]]             # keep relevant names
    names(par.wait) <- p                # use full names

    if (wait == "exponential")
    {
        if ("rate" %in% p)
            rate <- par.wait$rate
        else
            stop("parameter \"rate\" missing in 'par.wait'")
        m <- length(rate)

        if ("weights" %in% p && m > 1L)
            prob.w <- rep(par.wait$weights, length.out = m)
        else if (m == 1L)
            prob.w <- 1
        else
            stop("parameter \"weights\" missing in 'par.wait'")

        rates.w <- diag(-rate, m)
    }
    else if (wait == "Erlang")
    {
        if ("shape" %in% p)
            shape <- par.wait$shape
        else
            stop("parameter \"shape\" missing in 'par.wait'")
        if ("rate" %in% p)
            rate <- par.wait$rate
        else if ("scale" %in% p)
            rate <- 1 / par.wait$scale
        else
            stop("parameter \"rate\" or \"scale\" missing in 'par.wait'")
        if (length(shape) < length(rate))
            shape <- rep(shape, length.out = length(rate))
        else
            rate <-  rep(rate, length.out = length(shape))
        m <- sum(shape)

        if ("weights" %in% p && length(shape) > 1L)
        {
            prob.w <- numeric(sum(shape))
            prob.w[cumsum(c(1, head(shape, -1L)))] <- par.wait$weights
        }
        else if (length(shape) == 1L)
            prob.w <- c(1, rep(0, m - 1))
        else
            stop("parameter \"weights\" missing in 'par.wait'")

        rates.w <- diag(rep(-rate, shape), m)
        if (m > 1 && shape > 1)
        {
            tmp <- -head(diag(rates.w), -1L)
            tmp[cumsum(head(shape, -1L))] <- 0 # insert 0s in "ll corners"
            rates.w[cbind(seq_len(m - 1), seq(2, len = m - 1))] <- tmp
        }
    }
    else                                # wait == "phase-type"
    {
        if ("prob" %in% p)
            prob.w <- par.wait$prob
        else
            stop("parameter \"prob\" missing in 'par.wait'")
        if ("rates" %in% p)
            rates.w <- par.wait$rates
        else
            stop("parameter \"rates\" missing in 'par.wait'")
        m <- length(prob.w)
        if (!(is.matrix(rates.w) && nrow(rates.w) == m))
            stop("invalid parameters in 'par.wait'")
    }
    ## =================================================

    ## Empty definition of the output function. The body is set later.
    FUN <- function(u, survival = FALSE, lower.tail = !survival) {}

    ## Cramer-Lundberg model
    if (wait == "exponential" && m == 1L)
    {
        ## Special case with an explicit solution
        if (claims == "exponential" && n == 1L)
        {
            lambda <- -drop(rates.w) / premium.rate
            body(FUN) <- substitute({res <- a * exp(-(b) * u);
                                     if (lower.tail) res else 0.5 - res + 0.5},
                                    list(a = -lambda/drop(rates),
                                         b = -drop(rates) - lambda))
            environment(FUN) <- new.env() # new, empty environment
            class(FUN) <- c("ruin", class(FUN))
            return(FUN)
        }

        ## Use phase-type representation for all other claim severity models.
        pi <- drop(rates.w) * prob %*% solve(rates) / premium.rate
        Q <- rates - rowSums(rates) %*% pi
    }
    ## Sparre Andersen model (interarrival times other than single exponential)
    else
    {
        ## Matrix Q is a "fixed point" of some function (Asmussen &
        ## Rolski, 1992, p. 265-266). Many elements of this function
        ## never change, hence they are computed once and for all
        ## here.
        In <- diag(n)                    # n x n identity matrix
        Im <- diag(m)                    # m x m identity matrix
        t0pi <- -rowSums(rates) %o% prob # "multiple" of A(Q)
        A <- In %x% rbind(prob.w)        # first term of A(Q)
        B <- In %x% rates.w              # rhs of the Kronecker sum
        C <- In %x% -rowSums(rates.w)    # third term of A(Q)

        if (echo)
        {
            cat("Iteration\tMatrix Q (column major order)\n")
            exp <- expression(cat(" ", count, "\t\t ", Q1 <- Q,
                fill = TRUE))
        }
        else
            exp <- expression(Q1 <- Q)

        Q <- rates
        count <- 0L

        repeat
        {
            eval(exp)

            if (maxit < (count <- count + 1L))
            {
                warning("maximum number of iterations reached before obtaining convergence")
                break
            }

            Q1 <- Q
            Q <- rates - t0pi %*% A %*% solve(Q %x% Im + B, C)

            if (max(rowSums(abs(Q - Q1))) < tol)
                break
        }
        pi <- colSums(Q - rates) / (-sum(rates) * premium.rate)
    }

    ## Compute the probability of ruin using the cdf of a phase-type
    ## distribution with parameters pi and Q.
    body(FUN) <- substitute(pphtype(u, a, b, lower.tail = !lower.tail),
                            list(a = pi, b = Q))
    environment(FUN) <- new.env()       # new, empty environment
    class(FUN) <- c("ruin", class(FUN))
    FUN
}

plot.ruin <- function(x, from = NULL, to = NULL, add = FALSE,
                      xlab = "u", ylab = expression(psi(u)),
                      main = "Probability of Ruin", xlim = NULL, ...)
    curve(x, from = from, to = to, add = add, xlab = xlab,
          ylab = ylab, main = main, xlim = xlim, ...)
