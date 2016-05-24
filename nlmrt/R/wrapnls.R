wrapnls <- function(formula, start, trace = FALSE, data, lower = -Inf, 
    upper = Inf, control = list(), ...) {
    # A wrapper to call nlsmnq() and then call nls() with the
    # solution.  The calling sequence matches that of nlsmnq()
    # 
    # if ((is.null(lower) && is.null(upper))) { cat('NULL lower
    # and upper\n') lower<- (-Inf) upper<-Inf } if
    # ((all(lower)==(-Inf)) && (all(upper) == Inf)) cat('Inf
    # lower and upper\n') if (length(lower)==1) { cat('expand
    # lower\n') lower<-rep(lower, length(start)) } if
    # (length(upper)==1) { cat('expand upper\n')
    # upper<-rep(upper, length(start)) }
    
    if (is.null(data)) 
        stop("wrapnls() must have 'data' supplied")
    npar <- length(start)
    if (length(lower) < npar) {
        if (length(lower) == 1) 
            lower <- rep(lower, npar) else stop("lower bounds wrong length")
    }
    if (length(upper) < npar) {
        if (length(upper) == 1) 
            upper <- rep(upper, npar) else stop("upper bounds wrong length")
    }
    if (trace) {
        cat("wrapnls call with lower=")
        print(lower)
        cat("and upper=")
        print(upper)
    }
    # Note that there are no bounds or masks.
    first <- nlxb(formula, start, trace = trace, data = data, 
        lower = lower, upper = upper, control = control, ...)
    # Should check this has worked, but ...
    if (trace) 
        print(first)
    newstart <- first$coefficients
    names(newstart) <- names(start)
    # Should put this in a try(), but let's see it work first
    if (trace) {
        cat("newstart:")
        print(newstart)
    }
    if (all(is.infinite(lower)) && all(is.infinite(upper))) {
        if (trace) 
            cat("nls call with no bounds\n")
        second <- nls(formula, newstart, trace = trace, data = data, 
            control = control, ...)
    } else {
        if (trace) 
            cat("Now try nls - bounded\n")
        second <- nls(formula, newstart, trace = trace, data = data, 
            algorithm = "port", lower = lower, upper = upper, 
            control = control, ...)
    }
}
