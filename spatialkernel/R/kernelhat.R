## see notes at lambdahat()
## return NA at points outside polygon
kernelhat <- function(dpts, h, poly, win=apply(poly, 2, range),
                      delta=1, edge=TRUE)
{
    adapt <- chkernel()
    x <- seq(win[1,1]+delta/2, win[2,1], by=delta)
    y <- seq(win[1,2]+delta/2, win[2,2], by=delta)
    nx <- length(x)
    ny <- length(y)
    xygrids <- as.matrix(expand.grid(x=x, y=y))
    ndx <- which(pinpoly(poly, xygrids) > 0)
    ## 0 outside, 1 boundary, 2 inside, -1 error
    pts <- xygrids[ndx,]
    npts <- nrow(pts)
    ndpts <- nrow(dpts)
    res <- rep(NA, nx*ny)
    if(edge)
        c <- adaptpoly(pts, h, poly)$c
    else
        c <- rep(1, npts)
    ans <- .C("hat_lambda_c", as.double(pts), as.integer(npts),
              as.double(dpts), as.integer(ndpts), as.double(h),
              as.integer(adapt$kernel), as.double(c), lam=double(npts),
              PACKAGE="spatialkernel")$lam
    res[ndx] <- ans
    invisible(list(lambda=matrix(res, nrow=nx, ncol=ny), x=x, y=y,
                   dpts=dpts, h=h, poly=poly, win=win, delta=1, edge=edge))
}
