rghtsngv <- function(x, nv=min(nrow(x),ncol(x)), maxsqmatdim=10000)
{ 
	if (!is.matrix(x)) x <- as.matrix(x)
 	n <- nrow(x)
	p <- ncol(x)
	minnp <- min(n,p)	
	d <- numeric(minnp)
	vt <- numeric(minnp*p)
	retcode <- 999
	u <- numeric(n*p)
	iwork <- integer(8*minnp)
       	lwork <- max( 5*minnp+2*max(n,p),
		.Fortran("lworkgesdd",
			as.double(x),as.double(d),as.double(vt),
			as.integer(n),as.integer(p),as.integer(minnp),as.double(u),
			as.double(numeric(1)),as.integer(iwork),as.integer(retcode),
			PACKAGE="HiDimDA")[[8]]
		)
	work <- numeric(lwork)	
       	Fortout <- .Fortran("rsgvdgesdd",
			as.double(x),as.double(d),as.double(vt),
			as.integer(n),as.integer(p),as.integer(minnp),as.double(u),
			as.double(work),as.integer(lwork),as.integer(iwork),as.integer(retcode),
			PACKAGE="HiDimDA"
		)
	if (Fortout[[11]]!=0) {
		dgesddrtc <- Fortout[[11]]  
       		Fortout <- .Fortran("rsgvdgesvd",
			as.double(x),as.double(d),as.double(vt),
			as.integer(n),as.integer(p),as.integer(minnp),
			as.double(work),as.integer(lwork),as.integer(retcode),
			PACKAGE="HiDimDA"
		)
	   if (Fortout[[9]]!=0)  {
		dgesvdrtc <- Fortout[[9]]  
		if (p>maxsqmatdim) 
			stop(paste("Computation of right singular vectors failed for a",n,"by",p," matrix,\n  with return codes",dgesddrtc,"and",dgesvdrtc,"for the Lapack dgesdd and dgesvd functions.\n"))
		else {
			x2 <- numeric(p*p)
			d2 <- numeric(p)
       			Fortout <- .Fortran("rsgvdsyev",
				as.double(x),as.double(d),as.double(vt),as.double(x2),as.double(d2),
				as.integer(n),as.integer(p),as.integer(minnp),
				as.double(work),as.integer(lwork),as.integer(retcode),
				PACKAGE="HiDimDA"
			)
	   		if (Fortout[[11]]!=0) 
				stop(paste("Computation of right singular vectors failed for a",n,"by",p," matrix,\n  with return codes",dgesddrtc,dgesvdrtc,"and",Fortout[[11]],"for the Lapack dgesdd, dgesvd and dsyev functions.\n"))
		}
	   }
	}
	list(d=drop(Fortout[[2]]),u=NULL,v=t(matrix(Fortout[[3]],minnp,p))[,1:nv])
}

