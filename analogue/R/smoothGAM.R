## smoothGAM: smoother function supplied to prcurve
## wrapper to mcgv

`smoothGAM` <- function(lambda, x, choose = TRUE, complexity,
                        bs = "tp", ...,
                        family = gaussian(),
                        method = "REML",
                        select = FALSE,
                        control = list()) {
    ## complexity is the 'k' argument -
    ## choose selects whether to use fixed complexity or allow
    ## underlying fitting function to return complexity
    ord <- order(lambda)
    lambda <- lambda[ord]
    x <- x[ord]
    if(!missing(complexity)) {
        complexity <- round(complexity) ## move this out of smoothGAM
    } else {
        complexity <- -1
    }
    f <- gam(x ~ s(lambda, k = complexity, fx = choose, bs = bs),
             family = family, method = method, select = select,
             control = control, ...)
    p <- predict(f, x = lambda, type = "response")
    edf <- sum(f$edf[f$smooth[[1]]$first.para:f$smooth[[1]]$last.para]) + 1
    res <- list(lambda = lambda, x = x, fitted.values = p,
                complexity = edf, model = f)
    class(res) <- "prcurveSmoother"
    res
}
