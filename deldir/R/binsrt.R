binsrt <- function(x,y,rw) {
    n    <- length(x)
    ind  <- rep(0,n)
    rslt <- .Fortran("binsrt",
                     x=as.double(x),
                     y=as.double(y),
                     rw=as.double(rw),
                     n=as.integer(n),
                     ind=as.integer(ind),
                     rind=as.integer(ind),
                     tx=double(n),
                     ty=double(n),
                     ilst=integer(n),
                     nerror=as.integer(0),
                     PACKAGE="deldir"
                )
    if(rslt$nerror > 0)
        stop("Mismatch between number of points and number of sorted points.\n")
    list(x=rslt$tx,y=rslt$ty,ind=rslt$ind,rind=rslt$rind)
}
