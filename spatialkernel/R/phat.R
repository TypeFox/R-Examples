## {y}={0,1,2,..,m-1} converted inside
## calculate phat at point pts
phat<-function(gpts, pts, marks, h)
{
    adapt <- chkernel()
    ngpts <- length(gpts)/2 ## NO. of points
    n <- length(marks)
    ynames <- names(table(marks))
    m <- length(ynames)
    ynames0 <- 1:m -1 
    names(ynames0) <- ynames
    yy <- ynames0[as.character(marks)]
    c <- rep(1, ngpts)
    ans<-.C("hatpn", as.double(gpts), as.integer(ngpts), as.double(pts),
            as.integer(yy), as.integer(n), as.double(h), 
            as.integer(adapt$kernel),
            as.double(c), as.integer(m), p=double(ngpts*m),
            PACKAGE="spatialkernel")$p
    ans <- matrix(ans, ncol=m, dimnames=list(NULL, ynames))
    invisible(list(p=ans, pts=pts, gpts=gpts, marks=marks, h=h))
}
