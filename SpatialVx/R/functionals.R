ampstats <- function(X, Xhat, only.nonzero=FALSE) {
   if( only.nonzero) {
	idy <- Xhat != 0
	idx <- X != 0
   } else {
	xdim <- dim(Xhat)
	idy <- idx <- matrix( TRUE, xdim[1], xdim[2])
   }
   m1 <- mean( Xhat[idy], na.rm=TRUE)
   v1 <- var( c(Xhat[idy]), na.rm=TRUE)
   m2 <- mean( X[idx], na.rm=TRUE) 
   v2 <- var( c(X[idx]), na.rm=TRUE)
   if( !only.nonzero) {
	id <- !is.na( c(Xhat)) & !is.na( c(X))
	v12 <- cov( c(Xhat)[id], c(X)[id])
   } else v12 <- NA
   return( list( mean.fcst=m1, mean.vx=m2, var.fcst=v1, var.vx=v2, cov=v12))
} # end of 'ampstats' function.

UIQI <- function(X, Xhat, ...) {
   out <- list()
   data.name <- c(as.character(substitute(X)), as.character(substitute(Xhat)))
   names(data.name) <- c("verification","forecast")
   out$data.name <- data.name
   tmp <- ampstats( X=X, Xhat=Xhat, ...)
   sig12 <- sqrt( tmp$var.fcst*tmp$var.vx)
   if( !is.na( tmp$cov)) rho <- tmp$cov/sig12
   else rho <- NA
   bb <- 2*(tmp$mean.fcst*tmp$mean.vx)/(tmp$mean.fcst^2+tmp$mean.vx^2)
   dv <- 2*sig12/(tmp$var.fcst+tmp$var.vx)
   out$cor <- rho
   out$brightness.bias <- bb
   out$distortion.variability <- dv
   if( !is.na( rho)) out$UIQI <- rho*bb*dv
   else out$UIQI <- bb*dv
   return( out)
}

mae <- function(x1, x2, ...) {
   x1 <- c(x1)
   x2 <- c( x2)
   x <- matrix( abs( c(x1)-c(x2)), ncol=1)
   return( colMeans( x, na.rm=TRUE))
} # end of 'mae' function.

fft2d <- function(x, bigdim=NULL, ...) {
   bigdim.in <- !is.null( bigdim)
   if( !bigdim.in) out <- list()
   xdim <- dim( x)
   if( !bigdim.in) {
	bigdim <- xdim
	if( xdim[1] <= 1024) bigdim[1] <- 2^ceiling(log2(bigdim[1]))
 	else bigdim[1] <- ceiling( bigdim[1]/512)*512
	if( xdim[2] <= 1024) bigdim[2] <- 2^ceiling(log2(bigdim[2]))
	else bigdim[2] <- ceiling( bigdim[2]/512)*512
   } # end of if 'bigdim' not supplied stmts.
   hold <- matrix( 0, bigdim[1], bigdim[2])
   hold[1:xdim[1],1:xdim[2]] <- x
   hold <- fft( hold, ...)[1:xdim[1],1:xdim[2]]
   if( !bigdim.in) {
      out$fft <- hold
      out$bigdim <- bigdim
   } else out <- hold
   return(out)
} # end of 'fft2d' function.

surrogater2d <- function(Im, frac=0.95, n=10, lossfun="mae", maxiter=100, zero.down=TRUE, verbose=FALSE, ...) {
   xdim <- dim( Im)
   out <- array( NA, dim=c(xdim, n))
   Im.vec <- c( Im)
   tmp <- fft2d( Im)
   bigdim <- tmp$bigdim
   smN <- prod( xdim)
   bigN <- prod( bigdim)
   amps <- abs( tmp$fft)
   for( i in 1:n) {
      sur <- aaft2d( Im=Im, bigdim=bigdim)
   
      if( verbose) cat("\n", "Found initial surrogate via AAFT method.  Iterating to try to get same power spectrum and distribution as Im.\n")
      for( k in 1:maxiter) {
         if( verbose) cat(k, " ")
         temp <- fft2d( sur, bigdim=bigdim)
         ampsnew <- abs(temp)
         id <- ampsnew == 0
         if( any(id)) ampsnew[id] <- 1e-12
         ampsadj <- frac*amps + (1-frac)*ampsnew
         temp <- temp*ampsadj/ampsnew
	 # sur <- Re(fft2d( temp, bigdim=bigdim, inverse=TRUE)/bigN)
  	 sur <- Re(fft2d( temp, bigdim=bigdim, inverse=TRUE)/smN)
	 if( zero.down) {
	   sur[ sur < 0] <- 0
           sur <- zapsmall( sur)
	 } # end of if 'zero.down' stmt.
         rk0 <- rank( sur)
         S <- Im.vec[rk0]
         if(k == 1) {
   	   PastErr <- Inf
   	   Errmin <- Inf
         } # end of if k = 1 stmts.
         Err <- do.call( lossfun, c(list( x1=amps, x2=ampsnew), list(...)))
         if( Err < Errmin) {
   	if( verbose) cat("\n", "Err = ", Err, " < Errmin = ", Errmin, "\n")
           surmin <- sur
           Errmin <- Err
         } # end of if 'Err < Errmin' stmts.
         if( is.infinite( PastErr)) {
   	PastErr <- Err
   	next
         } else if( abs((Err-PastErr)/PastErr)*100 <= 1 ) break
	 PastErr <- Err
       } # end of for 'k' loop.
	if( verbose) cat("\n", "Finished iterating.  Final Err = ", Err, "\n")
	out[,,i] <- surmin
	if( verbose) cat("\n", "Found surrogate ", i, "\n")
    } # end of for 'i' loop.
    return( out)
} # end of 'surrogater2d' function.

aaft2d <- function(Im, bigdim=NULL) {
        xdim <- dim( Im)
        Im.vec <- c( Im)
        N <- length( Im.vec)
	Im2 <- matrix( rnorm( N), xdim[1], xdim[2])
        out <- Im2+NA
        Im2.vec <- c( Im2)
        rk1 <- rank( Im.vec)
        sIm2.vec <- sort( Im2.vec)
        Z <- sIm2.vec[rk1]
	if( is.null( bigdim)) {
           bigdim <- xdim
           if( xdim[1] <= 1024) bigdim[1] <- 2^ceiling(log2(bigdim[1]))
           else bigdim[1] <- ceiling( bigdim[1]/512)*512
           if( xdim[2] <= 1024) bigdim[2] <- 2^ceiling(log2(bigdim[2]))
           else bigdim[2] <- ceiling( bigdim[2]/512)*512
   	} # end of if 'bigdim' not supplied stmts.
        big0mat <- matrix( 0, bigdim[1], bigdim[2])
        bigN <- prod( bigdim)

        X <- big0mat
        X[1:xdim[1],1:xdim[2]] <- matrix(Z, xdim[1], xdim[2])
	Zfft <- fft2d(X, bigdim=bigdim)
        nr1 <- ceiling(xdim[1]/2)-1
        nr2 <- (xdim[1]-1)
        nc2 <- ceiling(xdim[2]/2)-1
        PI <- 2*pi*matrix( rnorm(nr1), ncol=1)
        PJ <- 2*pi*matrix( rnorm(nc2), nrow=1)
        PIJ <- 2*pi*matrix( rnorm(nr2*nc2), nr2, nc2)
        if( (xdim[2] %% 2) == 0) {
              QIJ <- 2*pi*matrix( rnorm(nr1), ncol=1)
              if( (xdim[1] %% 2) == 0) QIJ <- c( c(QIJ), 0, c(-QIJ[nr1:1,]))
              else QIJ <- c(QIJ, -QIJ[nr1:1,])
              PHIJ      <- cbind(PIJ, QIJ, -PIJ[nr2:1,nc2:1])
        } else PHIJ <- cbind( PIJ, -PIJ[nr2:1,nc2:1])
        if( (xdim[1] %% 2) == 0) PHI <- c(0, c(PI), 0, c(-PI[nr1:1,]))
        else PHI <- c(0, PI, -PI[nr1:1])
        if( (xdim[2] %% 2) == 0) PHJ <- c( 0, PJ, 0, -PJ[,nc2:1])
        else PHJ <- c( 0, c(PJ), c(-PJ[,nc2:1]))
        PH <- matrix(0, xdim[1], xdim[2])
        PH[2:xdim[1],2:xdim[2]] <- PHIJ
        PH[1,] <- PHJ
        PH[,1] <- PHI
        K <- big0mat
        K[1:xdim[1],1:xdim[2]] <- Zfft[1:xdim[1],1:xdim[2]]*exp(1i*PH)
        # J <- Re(fft(K, inverse=TRUE)/bigN)[1:xdim[1],1:xdim[2]]
	J <- Re(fft(K, inverse=TRUE)/N)[1:xdim[1],1:xdim[2]]
        rk2 <- rank( c( J))
        J <- sort( Im.vec)
        J <- J[ rk2]
        out <- matrix( J, xdim[1], xdim[2])
        return( out)
   } # end of 'aaft2d' function.

FQI <- function(object, surr=NULL, k=4, time.point=1, model=1, ...) {

    object <- locmeasures2dPrep(object=object, k=k)

    a <- attributes(object)
    out <- list()
    attributes(out) <- a

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(object, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(object, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(object, model=model)
    else dat <- datagrabber(object)

    X <- dat$X
    Y <- dat$Xhat

    dn <- a$data.name
    if(length(dn) == a$nforecast + 2) {
        mainname <- dn[1]
        dn <- dn[-1]
    } else mainname <- NULL

    vxname <- dn[1]
    dn <- dn[-1]

    if(!is.numeric(model)) {
        model.num <- (1:a$nforecast)[model == dn]
    } else model.num <- model
    attr(out, "data.name") <- c(mainname, vxname, dn[model.num])
    thresholds <- a$thresholds
    thresholds <- cbind(thresholds[,1], thresholds[,model.num + 1])
    colnames(thresholds) <- c(vxname, dn[model.num])
    attr(out, "thresholds") <- thresholds
    ## End: Get the data sets

    if( is.null( surr)) surr <- surrogater2d( Im=X, ...)
    xdim <- dim( X)
    Yim <- im( Y)
    Xim <- im( X)
    q <- dim( thresholds)[1]
    ks <- a$k
    nk <- length( ks)
    phd.norm <- fqi <- matrix( NA, nrow=nk, ncol=q)
    uiqi.norm <- numeric( q)+NA
    locperfer <- function(x1, x2, thresh, k, ...) {
 	x1 <- im( x1)
 	Ix1 <- solutionset( x1 >= thresh)
	return( locperf( X=x2, Y=Ix1, which.stats="ph", k=k, ...)$ph)
    } # end of 'locperfer' internal function.
    for( threshold in 1:q) {
	Ix <- solutionset( Xim >= thresholds[ threshold,1])
	Iy <- solutionset( Yim >= thresholds[ threshold,2])
	idx <- X >= thresholds[ threshold,1]
	idy <- Y >= thresholds[ threshold,2]
	X2 <- X
	X2[!idx] <- 0
	Y2 <- Y
	Y2[!idy] <- 0
	denom <- UIQI( X=X2, Xhat=Y2, only.nonzero=TRUE)$UIQI
	uiqi.norm[threshold] <- denom
	for( k in 1:nk) {
	   num1 <- locperf( X=Ix, Y=Iy, which.stats="ph", k=ks[k])$ph
	   num2 <- mean( apply(surr, 3, locperfer, x2=Ix, thresh=thresholds[threshold,1], k=ks[k]), na.rm=TRUE)
	   num.tmp <- num1/num2
	   phd.norm[k,threshold] <- num.tmp
	   fqi[k,threshold] <- num.tmp/denom
	} # end of for 'k' loop.
    } # end of for 'threshold' loop.
    out$phd.norm <- phd.norm
    out$uiqi.norm <- uiqi.norm
    out$fqi <- fqi
    class(out) <- "fqi"
    return(out)
} # end of 'FQI' function.

summary.fqi <- function(object, ...) {
   x <- attributes(object)
   if( !is.null( x$qs)) lu <- x$qs
   else {
	u <- x$thresholds
	q <- dim( u)[1]
	if( all( u[,1] == u[,2])) lu <- as.character(u[,1])
	else lu <- paste(1:q, ": (", u[,1], ", ", u[,2], "); ", sep="")
   }
   cat("Comparison for:\n ", x$data.name, "\n")
   cat("\n", "Thresholds are: ", lu, "\n")
   k <- paste("k = ", as.character( x$k), "; ", sep="")

   y <- object$phd.norm
   rownames( y) <- k
   colnames( y) <- lu
   cat("\n", "normalized Partial Hausdorff distance\n")
   print( y)

   y <- matrix( object$uiqi.norm, nrow=1)
   colnames( y) <- lu
   cat("\n", "modified universal image quality index\n")
   print( y)

   y <- object$fqi
   rownames( y) <- k
   colnames( y) <- lu
   cat("\n", "Forecast Quality Index (FQI)\n")
   print( y)

   invisible()
} # end of 'summary.fqi' function.

print.fqi <- function(x, ...) {

    print(summary(x))
    invisible()
} # end of 'print.fqi' function.
