LORgee.control <-
function (tolerance = 0.001, maxiter = 15, verbose = FALSE, TRACE = FALSE) 
{
    if (!is.numeric(tolerance) || tolerance <= 0) 
        stop("value of LORgee's 'tolerance' must be > 0")
    if (!is.numeric(maxiter) || maxiter <= 0) 
        stop("maximum number of LORgee's iterations must be > 0")
    list(tolerance = tolerance, maxiter = maxiter, verbose = verbose, 
        TRACE = TRACE)
}

