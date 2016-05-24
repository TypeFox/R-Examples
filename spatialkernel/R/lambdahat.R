## kernel density estimation of \hat\lambda
## should be divided by |A| so that N(A) follows a Poisson with mean lambda*|A|
## See Peter (2003), pp47 homogeneous Poisson process definition
lambdahat<-function(pts, h, gpts=NULL, poly=NULL, edge=TRUE)
{
  adapt <- chkernel()
  npts <- nrow(pts)
  if(is.null(gpts)) { ##Baddeley's modified version
	if(edge) c <- adaptpoly(pts, h, poly)$c else c <- rep(1, npts)
    ans <- .C("hat_lambda_b", as.double(pts), as.integer(npts),
            as.double(h), as.integer(adapt$kernel), as.double(c),
            lam=double(npts), PACKAGE="spatialkernel")$lam
  } else {
    ngpts <- nrow(gpts)
    if(edge) c <- adaptpoly(gpts, h, poly)$c else c <- rep(1, ngpts)
    ans <- .C("hat_lambda_c", as.double(gpts), as.integer(ngpts),
            as.double(pts), as.integer(npts), as.double(h),
            as.integer(adapt$kernel), as.double(c), lam=double(ngpts),
            PACKAGE="spatialkernel")$lam
  }
  invisible(list(lambda=ans, pts=pts, gpts=gpts, poly=poly, h=h, edge=edge))
}

adaptpoly<-function(pts, h, poly) 
{
    c1 <- 10; eps <- 0.0001; mcalls <- 10000
    adapt <- chkernel()
	if(adapt$kernel==1) c2 <- 20 else c2 <- 1
    rng <- c(range(poly[,1]), range(poly[,2]))
    if(is.null(nrow(pts))) npts <- 1 else npts <- nrow(pts)
    ans<-.C("adaptpoly", as.double(poly), as.integer(nrow(poly)), as.double(pts),
            as.integer(npts), as.double(h), as.integer(adapt$kernel),
            as.double(c1), as.double(c2), as.double(rng),
            as.double(eps), err=double(npts), as.integer(mcalls),
            ncalls=integer(npts), ier=integer(6), c=double(npts),
			PACKAGE="spatialkernel")
    invisible(list(c=ans$c, err=ans$err, ncalls=ans$ncalls, ier=ans$ier,
        pts=pts, h=h, poly=poly))
}
