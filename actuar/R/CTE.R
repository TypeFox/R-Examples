### ===== actuar: An R Package for Actuarial Science =====
###
### Conditional Tail Expectation for objects of class 'aggregateDist'.
###
### AUTHORS: Tommy Ouellet, Vincent Goulet <vincent.goulet@act.ulaval.ca>

CTE <- function(x, ...)
    UseMethod("CTE")

CTE.aggregateDist <- function(x, conf.level = c(0.9, 0.95, 0.99),
                              names = TRUE, ...)
{
    chkDots(...)                        # method does not use '...'

    label <- comment(x)

    ## Normal approximation; an exact formula is available
    if (label == "Normal approximation")
    {
        m <- get("mean", environment(x))
        sd <- sqrt(get("variance", environment(x)))
        res <- m + sd * dnorm(qnorm(conf.level)) / (1 - conf.level)
    }
    ## Normal Power approximation; explicit formula in Castaner,
    ## Claramunt and Marmol (2013)
    else if (label == "Normal Power approximation")
    {
        m <- get("mean", envir = environment(x))
        sd <- sqrt(get("variance", envir = environment(x)))
        sk <- get("skewness", envir = environment(x))

        q <- qnorm(conf.level)
        res <- m + sd * dnorm(q) * (1 + sk * q/6) / (1 - conf.level)
    }
    ## Recursive method, simulation and convolutions; each yield a
    ## step function that can be used to make calculations.
    else
    {
        val <- get("x", envir = environment(x))
        prob <- get("fs", envir = environment(x))
        f2 <- function(a)
        {
            pos <- val > VaR(x, a)
            drop(crossprod(val[pos], prob[pos])) / sum(prob[pos])
        }
        res <- sapply(conf.level, f2)
    }

    if (names)
    {
        dig <- max(2, getOption("digits"))
        names(res) <- formatC(paste(100 * conf.level, "%", sep = ""),
                              format = "fg", width = 1, digits = dig)
    }
    res
}

TVaR <- CTE
