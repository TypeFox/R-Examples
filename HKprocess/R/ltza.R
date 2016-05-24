ltza <- function(r,x) {
    nx <- as.integer(length(x))
    EPS <- .Machine$double.eps
    # Call the ltza.c from the C library

    out <- .C("ltza",as.double(r),nx,as.double(x),nx,EPS,
    tr = array(0,dim = c(1,4)),fault = as.integer(1),PACKAGE = "HKprocess")

    return(setNames(as.vector(out$tr),c("t(x) * solve(R) * x",
    "t(en) * solve(R) * x","t(en) * solve(R) * en",
    "natural logarithm of determinant")))
}