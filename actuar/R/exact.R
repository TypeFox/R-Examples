### ===== actuar: An R Package for Actuarial Science =====
###
### Exact calculation of the aggregate claim amount distribution
### function by convolution. Requires a discrete distribution for
### claim amounts.
###
### AUTHORS:  Vincent Goulet <vincent.goulet@act.ulaval.ca>
### and Louis-Philippe Pouliot

exact <- function(fx, pn, x.scale = 1)
{
    ## Some useful lengths
    m <- length(fx)                   # 1 + maximum claim amount
    n <- length(pn) - 1               # maximum number of claims
    r <- n * m - n + 1                # maximum total amount of claims

    ## Initialization of the output vector
    fs <- rep(0, r)
    fs[1] <- pn[1]                    # Pr[N = 0]

    ## Convolutions
    fxc <- 1
    for (i in 1:n)
    {
        pos <- seq_len(i * m - i + 1)
        fxc <- convolve(fx, rev(fxc), type = "open")
        fs[pos] <- fs[pos] + fxc * pn[i + 1]
    }

    FUN <- approxfun((0:(length(fs) - 1)) * x.scale, pmin(cumsum(fs), 1),
                     method = "constant", yleft = 0, yright = 1, f = 0,
                     ties = "ordered")
    class(FUN) <- c("ecdf", "stepfun", class(FUN))
    assign("fs", fs, envir = environment(FUN))
    assign("x.scale", x.scale, envir = environment(FUN))
    FUN
}
