ltzb <- function(r,x) {
    nx <- as.integer(length(x))
    EPS <- .Machine$double.eps

    out <- .C("ltzb",as.double(r),nx,as.double(x),nx,EPS,
    tr = array(0,dim=c(1,2)),fault = as.integer(1),PACKAGE = "HKprocess")

    return(setNames(as.vector(out$tr),c("t(en) * solve(R) * x",
    "t(en) * solve(R) * en")))
}