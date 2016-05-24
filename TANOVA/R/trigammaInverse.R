trigammaInverse<-function (x) 
{
    if (!is.numeric(x)) 
	stop("Non-numeric argument to mathematical function")
    if (length(x) == 0) 
	return(numeric(0))
    omit <- is.na(x)
    if (any(omit)) {
        y <- x
        if (any(!omit)) 
		y[!omit] <- Recall(x[!omit])
        return(y)
    }
    omit <- (x < 0)
    if (any(omit)) {
        y <- x
        y[omit] <- NaN
        warning("NaNs produced")
        if (any(!omit)) 
		y[!omit] <- Recall(x[!omit])
        return(y)
    }
    omit <- (x > 1e+07)
    if (any(omit)) {
        y <- x
        y[omit] <- 1/sqrt(x[omit])
        if (any(!omit)) 
		y[!omit] <- Recall(x[!omit])
        return(y)
    }
    omit <- (x < 1e-06)
    if (any(omit)) {
        y <- x
        y[omit] <- 1/x[omit]
        if (any(!omit)) 
		y[!omit] <- Recall(x[!omit])
        return(y)
    }
    y <- 0.5 + 1/x
    iter <- 0
    repeat {
        iter <- iter + 1
        tri <- trigamma(y)
        dif <- tri * (1 - tri/x)/psigamma(y, deriv = 2)
        y <- y + dif
        if (max(-dif/y) < 1e-08) 
		break
        if (iter > 50) {
            warning("Iteration limit exceeded")
            break
        }
    }
    y
}
