GetFitHReg <-
function (y, t = 1:length(y),nf=150) 
{
    if (length(y) < 6)
        stop("length of y must be 6 or more")
    if (!(is.numeric(y)&&is.numeric(t)))
        stop("y and t must be numeric")
    if (length(y)!=length(t))
        stop("length(y)!=length(t)")
    theta <- numeric(2)
    theta[1] <- length(y)
    theta[2] <- nf #Number of frequencies enumerated
    ans <- .C("GetHReg", y = as.double(y), t = as.double(t), 
        theta = as.double(theta), PACKAGE = "pRSR")
    ans$theta
}
