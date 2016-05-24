### ===== actuar: An R Package for Actuarial Science =====
###
### Simulation of a aggregate claim amounts
###
### AUTHORS:  Vincent Goulet <vincent.goulet@act.ulaval.ca>
### and Louis-Philippe Pouliot

simS <- function(n, model.freq, model.sev)
{
    ## Prepare the call to simul() by building up 'nodes'
    level.names <- names(if (is.null(model.freq)) model.sev else model.freq)
    nlevels <- length(level.names)
    nodes <- as.list(c(rep(1, nlevels - 1), n))
    names(nodes) <- level.names

    ## Get sample
    x <- aggregate(simul(nodes = nodes,
                         model.freq = model.freq,
                         model.sev = model.sev))[-1]

    ## Compute the empirical cdf of the sample. Done manually instead
    ## of calling stats:::ecdf() to keep a copy of the empirical pmf
    ## in the environment without computing it twice.
    x <- sort(x)
    vals <- unique(x)
    fs <- tabulate(match(x, vals))/length(x)
    FUN <- approxfun(vals, pmin(cumsum(fs), 1), method = "constant",
                     yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(FUN) <- c("ecdf", "stepfun", class(FUN))
    assign("fs", fs, envir = environment(FUN))
    FUN
}
