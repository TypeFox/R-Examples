## smoothSpline: smootherFcn supplied to prcurve
## again a wrapper but allow us to specify De'ath's recommended
## smoother strategy
`smoothSpline` <- function(lambda, x, choose = TRUE,
                           complexity, ..., penalty = 1,
                           cv = FALSE, keep.data = FALSE,
                           control.spar = list(low = 0)) {
    ## complexity is the 'df' argument
    ## choose selects whether to use fixed complexity or allow
    ## underlying fitting function to return complexity
    ord <- order(lambda)
    lambda <- lambda[ord]
    x <- x[ord]
    if(choose) { ## choose complexity
        f <- smooth.spline(lambda, x, ...,
                           penalty = penalty,
                           keep.data = keep.data, cv = cv,
                           control.spar = control.spar)
    } else { ## use specified complexity
        f <- smooth.spline(lambda, x, ..., df = complexity,
                           penalty = penalty, ## no cv as specifying df
                           keep.data = keep.data,
                           control.spar = control.spar)
    }
    p <- predict(f, x=lambda)$y
    res <- list(lambda = lambda, x = x, fitted.values = p,
                complexity = f$df, model = f)
    class(res) <- "prcurveSmoother"
    res
}
