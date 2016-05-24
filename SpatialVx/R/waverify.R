makeEmptyWave <- function(x) {
   n <- length(x)
   for(i in 1:n) x[[i]] <- x[[i]]*0
   return(x)
} # end of 'makeEmptyWave' function.

makeWaveNames <- function(J) {
   out <- list()
   out$HH <- paste( "HH", 1:J, sep="")
   out$HL <- paste( "HL", 1:J, sep="")
   out$LH <- paste( "LH", 1:J, sep="")
   out$LL <- paste( "LL", J, sep="")
   return(out)
} # end of 'makeWaveNames' function.

waverify2d <- function(X, ..., Clim=NULL, wavelet.type="haar", J=NULL) {
    UseMethod("waverify2d", X)
} # end of 'waverify2d' function.

waverify2d.SpatialVx <- function(X, ..., Clim=NULL, wavelet.type="haar", J=NULL,
    useLL=FALSE, compute.shannon=FALSE, which.space="field", time.point=1, model=1, verbose=FALSE) {

    out <- list()
    a <- attributes(X)
    attributes(out) <- a

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(X, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(X, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(X, model=model)
    else dat <- datagrabber(X)

    x <- dat$X
    y <- dat$Xhat
    ## End: Get the data sets

    out2 <- waverify2d.default(X=x, Y=y, Clim=Clim, wavelet.type=wavelet.type, J=J,
			        useLL=useLL, compute.shannon=compute.shannon,
				which.space=which.space, verbose=verbose)

    out$A <- out2$A
    out$J <- out2$J
    out$X.wave <- out2$X.wave
    out$Y.wave <- out2$Y.wave
    out$Clim.wave <- out2$Clim.wave
    out$Shannon.entropy <- out2$Shannon.entropy
    out$energy <- out2$energy
    out$mse <- out2$mse
    out$rmse <- out2$rmse
    out$acc <- out2$acc
    out$wave.method <- out2$wave.method

    attr(out, "time.point") <- time.point
    attr(out, "model") <- model
    class(out) <- "waverify2d"
    return(out)
} # end of 'waverify2d.SpatialVx' function.

waverify2d.default <- function(X, ..., Y, Clim=NULL,
			wavelet.type="haar", J=NULL, useLL=FALSE,
			compute.shannon=FALSE, which.space="field",
			verbose=FALSE) {
   ##
   ## Function to apply wavelet decomposition to each field of a verification
   ## set, V=(Y, X), containing the forecast and observation/analysis fields,
   ## respectively.  Various verification metrics are applied to each wave number.
   ##
   ## 'X', 'Y', 'Clim'  an observed and forecast field, and optional climate field (to calculate the ACC).
   ## 'wavelet.type' character name giving the wavelet decomposition function to be used
   ##	(see help file for 'dwt.2d' from 'waveslim' package).
   ## 'J' numeric giving the number of levels of the wavelet to use.  If NULL, uses log2(dim( Y)[1]).
   ## 'useLL' logical, if TRUE the average component "LL" is used for the inverse wavelet transform, which ensures positive energy.
   ## 'compute.shannon' logical whether or not to compute the Shannon entropy.
   ## 'verbose' logical, should progress information be printed to the screen?
   ##
   ## Details: A diadic wavelet (in future, may also allow nondiadic fields) decomposition is applied to both forecast and observed
   ##	fields.  The rmse is taken for each scale of the decomposition.  If a climatology field is also supplied,
   ##	then the ACC is also take for each scale.
   ##
   ## Value: a list object containing the components:
   ##
   ## 'Y.wave', 'X.wave' objects of class "dwt.2d" or "modwt.2d" (if fields are not dyadic) giving the
   ##	wavelet decomposition coefficients for each field.
   ## 'J' the number of levels used for the wavelet decomposition.
   ## 'Shannon.entropy' 2X1 matrix (if compute.shannon is TRUE) giving the Shannon entropy score for the chosen wavelet function applied to each image.
   ## 
   ## Depends: pacakges 'waveslim' 
   ##
   if( verbose)  begin.time <- Sys.time()
   out <- list()
    attributes(out) <- atmp <- list(...)
    if(is.null(atmp$loc.byrow)) attr(out, "loc.byrow") <- FALSE

   climate <- !is.null( Clim)
   
   N <- dim( Y)[1]
   
   ## Do the 2-d wavelet transforms.
   if( is.null( J)) J <- log2(N)
   out$J <- J
   
   wv.X <- dwt.2d( X, wf=wavelet.type, J=J)
   wv.Y <- dwt.2d( Y, wf=wavelet.type, J=J)
   if( climate) wv.Clim <- dwt.2d( Clim, wf=wavelet.type, J=J)
   out$Y.wave <- wv.Y
   out$X.wave <- wv.X
   if( climate) out$Clim.wave <- wv.Clim
   
   if(compute.shannon) {
      Shannon.entropy <- matrix( NA, nrow=1, ncol=2)
      colnames( Shannon.entropy) <- c("Verification", "Forecast")
      tmp <- abs(unlist(wv.X))
      tmp <- tmp/sum(tmp,na.rm=TRUE)
      id0 <- tmp==0
      tmp[!id0] <- tmp[!id0]*log(tmp[!id0])
      Shannon.entropy[,1] <- -sum(tmp,na.rm=TRUE)
      tmp <- abs(unlist(wv.Y))
      tmp <- tmp/sum(tmp,na.rm=TRUE)
      id0 <- tmp==0
      tmp[!id0] <- tmp[!id0]*log(tmp[!id0])
      Shannon.entropy[,2] <- -sum(tmp, na.rm=TRUE)
      out$Shannon.entropy <- Shannon.entropy
   }
   
   # List names for coefficients of each wavelet object.
   tmp <- makeWaveNames(J=J)
   HH <- tmp$HH
   HL <- tmp$HL
   LH <- tmp$LH
   LL <- tmp$LL

   mse.i <- rmse.i <- numeric(J) + NA
   if( climate) acc.i <- rmse.i
   
   if( !climate) {
      energy.i <- matrix( NA, nrow=J, ncol=2)
      colnames( energy.i) <- c("Verification", "Forecast")
   } else {
      energy.i <- matrix( NA, nrow=J, ncol=3)
      colnames( energy.i) <- c("Verification", "Forecast", "Climatology")
   }
  
   if(verbose) cat("\n", "calculating statistics at each level.\n")
   for( i in 1:J) {
     if( verbose) cat(i, " ")
     energy.i[i,1] <- sum(wv.X[[HH[i]]]^2, na.rm=TRUE) + sum(wv.X[[HL[i]]]^2, na.rm=TRUE) + sum(wv.X[[LH[i]]]^2, na.rm=TRUE)
     energy.i[i,2] <- sum(wv.Y[[HH[i]]]^2, na.rm=TRUE) + sum(wv.Y[[HL[i]]]^2, na.rm=TRUE) + sum(wv.Y[[LH[i]]]^2, na.rm=TRUE)
     if( climate) energy.i[i,3] <- sum(wv.Clim[[HH[i]]]^2, na.rm=TRUE) + sum(wv.Clim[[HL[i]]]^2, na.rm=TRUE) + sum(wv.Clim[[LH[i]]]^2, na.rm=TRUE)
   
      if(verbose) {
	   if(which.space=="field") cat("Finding inverse transforms for this scale (i.e., detail reconstruction).\n")
	   else cat("Finding coefficient detail for this scale.\n")
      }
      inv.hold0.i <- detailer(wv.X, level=i, which.space=which.space, useLL=useLL, lnames=tmp)
      inv.hold1.i <- detailer(wv.Y, level=i, which.space=which.space, useLL=useLL, lnames=tmp)
      if( climate) inv.hold2.i <- detailer(wv.Clim, level=i, which.space=which.space, useLL=useLL, lnames=tmp)
      if(verbose) cat("Inverse transforms for present scale found.  Calculating vx stats.\n")
      mse.i[i] <- mean( (c( inv.hold0.i) - c( inv.hold1.i))^2, na.rm=TRUE)
      rmse.i[i] <- sqrt( mse.i[i])
   
      if( climate) {
   	if(verbose) cat("Calculating ACC.\n")
   	denom <- sqrt( sum( (c( inv.hold0.i) - c( inv.hold2.i))^2, na.rm=TRUE))*sqrt( sum( (c( inv.hold1.i) - c( inv.hold2.i))^2, na.rm=TRUE))
   	numer <- sum( diag( t( inv.hold0.i - inv.hold2.i) %*% (inv.hold1.i - inv.hold2.i)), na.rm=TRUE)
   	acc.i[i] <- numer/denom
   	if(verbose) cat("ACC calculated.\n")
      }
   } # end of for 'i' loop.
   
   out$energy <- energy.i
   out$mse <- mse.i
   out$rmse <- rmse.i
   if( climate) out$acc <- acc.i
   out$type <- wavelet.type
   out$wave.method <- "DWT"
   if(verbose) print( Sys.time() - begin.time)
   class(out) <- "waverify2d"
   return(out)
} # end of 'waverify2d.default' function.

print.waverify2d <- function(x, ...) {
    a <- attributes(x)
    if(!is.null(a$data.name)) {
	cat(a$data.name, "\n")
	cat(a$msg)
    }

    cat("\n", "Wavelet method = ", x$wave.method)
    cat("\n", " number of levels (J) = ", x$J, "\n")
    cat("Energy:\n")
    print(x$energy)
    cat("\n", "MSE:\n")
    print(x$mse)
    cat("\n", "RMSE:\n")
    print(x$rmse)

    invisible()
} # end of 'print.waverify2d' function.

plot.waverify2d <- function(x, ..., main1="X", main2="Y", main3="Climate",
    which.plots=c("all", "dwt2d", "details", "energy", "mse", "rmse", "acc"),
    separate=FALSE, col, horizontal=TRUE) {

    if(missing(col)) col <- c("gray", tim.colors(64))

    which.plots <- tolower(which.plots)
    which.plots <- match.arg(which.plots)

    a <- attributes(x)
    loc.byrow <- a$loc.byrow

    if(!is.null(a$projection) && any(is.element(c("all","details"), which.plots))) {
        proj <- a$projection
	domap <- a$map
	loc <- a$loc
	xd <- a$xdim
	if(proj) loc2 <- list(x=matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow), y=matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow))
    } else proj <- domap <- FALSE

    if(domap) {
	locr <- apply(loc, 2, range, finite=TRUE)
	ax <- list(x=pretty(round(loc[,1], digits=2)), y=pretty(round(loc[,2], digits=2)))
    }

    if(!is.null(a$data.name)) {
	if(length(a$data.name) == a$nforecast + 2) dn <- a$data.name[-1]
	else dn <- a$data.name
	if(is.numeric(a$model)) model.num <- a$model
	else model.num <- (1:a$nforecast)[a$model == dn[-1]]
	dn <- c(dn[1], dn[model.num + 1])
    }

    if(missing(main1) && !is.null(a$data.name)) main1 <- dn[1]
    if(missing(main2) && !is.null(a$data.name)) main2 <- dn[2]


    if(!is.null(a$msg)) par(oma=c(0,0,2,0))

    X <- x$X.wave
    Y <- x$Y.wave
    if(is.null(X) | is.null(Y)) do.detail <- FALSE
    else do.detail <- TRUE
    if(climate <- !is.null(x$Clim.wave)) Clim <- x$Clim.wave

    if(any(is.element(c("all","dwt2d2"),which.plots)) & (x$wave.method == "DWT")) {
	if(!do.detail) stop("plot.waverify2d: sorry, no detail information provided for one or both of the fields.")
	if(!separate) {
	   if(climate) nc <- 3
	   else nc <- 2
	   par(mfrow=c(1,nc), mar=c(4,4,4,4))
	}
	plot(X, main=main1, ...)
	plot(Y, main=main2, ...)
	if(climate) plot(Clim, main=main3, ...)
    }
    J <- x$J
    if(any(is.element(c("all", "details"), which.plots))) {
      if(!do.detail) stop("plot.waverify2d: sorry, no detail information provided for one or both of the fields.")
      tmp <- makeWaveNames(J=J)
      HH <- tmp$HH
      HL <- tmp$HL
      LH <- tmp$LH
      LL <- tmp$LL
      if(!separate) {
         nr <- ceiling(J/2)
         if(climate) nc <- 6
         else nc <- 4
         par(mfrow=c(nr,nc), mar=c(1,2,4,1))
      }
      holdX <- makeEmptyWave(X)
      holdY <- holdX
      if(climate) holdClim <- holdY
      for(i in 1:J) {
	ni <- length(c(holdX[[HH[i]]]))*3
	holdX[[HH[i]]] <- X[[HH[i]]]
	holdX[[HL[i]]] <- X[[HL[i]]]
	holdX[[LH[i]]] <- X[[LH[i]]]
	holdY[[HH[i]]] <- Y[[HH[i]]]
        holdY[[HL[i]]] <- Y[[HL[i]]]
        holdY[[LH[i]]] <- Y[[LH[i]]]
	if(climate) {
	   holdClim[[HH[i]]] <- Clim[[HH[i]]]
           holdClim[[HL[i]]] <- Clim[[HL[i]]]
           holdClim[[LH[i]]] <- Clim[[LH[i]]]
	}

	lamX <- sqrt( sum( c( holdX[[HH[i]]]^2) + c(holdX[[HL[i]]]^2) + c(holdX[[LH[i]]]^2), na.rm=TRUE)/ni)
	holdX[[HH[i]]] <- sign(X[[HH[i]]])*pmax(c(abs(X[[HH[i]]]) - lamX), 0)
        holdX[[HL[i]]] <- sign(X[[HL[i]]])*pmax(c(abs(X[[HL[i]]]) - lamX), 0)
        holdX[[LH[i]]] <- sign(X[[LH[i]]])*pmax(c(abs(X[[LH[i]]]) - lamX), 0)

	lamY <- sqrt(sum(c(holdY[[HH[i]]]^2) + c(holdY[[HL[i]]]^2) + c(holdY[[LH[i]]]^2), na.rm=TRUE)/ni)
        holdY[[HH[i]]] <- sign(Y[[HH[i]]])*pmax(c(abs(Y[[HH[i]]]) - lamY), 0)
        holdY[[HL[i]]] <- sign(Y[[HL[i]]])*pmax(c(abs(Y[[HL[i]]]) - lamY), 0)
        holdY[[LH[i]]] <- sign(Y[[LH[i]]])*pmax(c(abs(Y[[LH[i]]]) - lamY), 0)

	if(climate) {
	   lamClim <- sqrt( sum( c( holdClim[[HH[i]]]^2) + c(holdClim[[HL[i]]]^2) + c(holdClim[[LH[i]]]^2), na.rm=TRUE)/ni)
           holdClim[[HH[i]]] <- sign(Clim[[HH[i]]])*pmax(c(abs(Clim[[HH[i]]]) - lamClim), 0)
           holdClim[[HL[i]]] <- sign(Clim[[HL[i]]])*pmax(c(abs(Clim[[HL[i]]]) - lamClim), 0)
           holdClim[[LH[i]]] <- sign(Clim[[LH[i]]])*pmax(c(abs(Clim[[LH[i]]]) - lamClim), 0)
	}
	if(x$wave.method=="DWT") {
	   invX <- idwt.2d(holdX)
	   invY <- idwt.2d(holdY)
	   if(climate) invClim <- idwt.2d(holdClim)
	} else {
	   invX <- imodwt.2d(holdX)
	   invY <- imodwt.2d(holdY)
	   if(climate) invClim <- imodwt.2d(holdClim)
	}

	if(!is.null(x$loc)) loc <- x$loc
	else {
	   xdim <- dim(invX)
	   loc <- cbind(rep(1:xdim[1],xdim[2]),rep(1:xdim[2],each=xdim[1]))
	}
	if(!climate) zl <- range(c(c(invX),c(invY)),na.rm=TRUE)
	else zl <- range(c(c(invX),c(invY),c(invClim)),na.rm=TRUE)

	if(domap) {
	    if(proj) {
		map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
		poly.image(loc2$x, loc2$y, invX, main=paste(main1, " (level ", i, ")", sep=""),
				zlim=zl, col=col, axes=FALSE, add=TRUE, ...)
		map(add=TRUE, lwd=2)
		map(add=TRUE, database="state")

		map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
		poly.image(loc2$x, loc2$y, invY, main=paste(main2, " (level ", i, ")", sep=""), zlim=zl, col=col, axes=FALSE, add=TRUE, ...)
		map(add=TRUE, lwd=2)
                map(add=TRUE, database="state")

		if(climate) {
                    map(xlim=locr[,1], ylim=locr[,2], type="n")
		    axis(1, at=ax$x, labels=ax$x)
		    axis(2, at=ax$y, labels=ax$y)
		    poly.image(loc2$x, loc2$y, invClim, main=paste(main3, " (level ", i, ")", sep=""), zlim=zl, axes=FALSE, add=TRUE, ...)
		    map(add=TRUE, lwd=2)
                    map(add=TRUE, database="state")
		}
		title("")
	    } else {
		map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
		image(as.image(invX, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main=paste(main1, " (level ", i, ")", sep=""),
                                zlim=zl, col=col, axes=FALSE, add=TRUE, ...)
		map(add=TRUE, lwd=2)
                map(add=TRUE, database="state")

		map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
		image(as.image(invY, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main=paste(main2, " (level ", i, ")", sep=""),
				zlim=zl, col=col, axes=FALSE, add=TRUE, ...)
		map(add=TRUE, lwd=2)
                map(add=TRUE, database="state")

		if(climate) {
		    map(xlim=locr[,1], ylim=locr[,2], type="n")
		    axis(1, at=ax$x, labels=ax$x)
		    axis(2, at=ax$y, labels=ax$y)
		    image(as.image(invClim, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE),
				main=paste(main3, " (level ", i, ")", sep=""), zlim=zl, axes=FALSE, add=TRUE, ...)
		    map(add=TRUE, lwd=2)
                    map(add=TRUE, database="state")
		}
		title("")
	    }

	} else {

	    if(proj) {
		poly.image(loc2$x, loc2$y, invX, main=paste(main1, " (level ", i, ")", sep=""), zlim=zl, col=col, axes=FALSE, ...)
		poly.image(loc2$x, loc2$y, invY, main=paste(main2, " (level ", i, ")", sep=""), zlim=zl, col=col, axes=FALSE, ...)
		if(climate) poly.image(loc2$x, loc2$y, invClim, main=paste(main3, " (level ", i, ")", sep=""), zlim=zl, col=col, axes=FALSE, ...)
	    } else {
	        image(invX, main=paste(main1, " (level ", i, ")", sep=""), zlim=zl, col=col, axes=FALSE, ...)
	        image(invY, main=paste(main2, " (level ", i, ")", sep=""), zlim=zl, col=col, axes=FALSE, ...)
	        if(climate) image(invX, main=paste(main3, " (level ", i, ")", sep=""), zlim=zl, col=col, axes=FALSE, ...) 
	    }
	    title("")
	} # end of if else 'domap' stmts.
	image.plot(invX, legend.only=TRUE, horizontal=horizontal, col=col, zlim=zl, ...)

	holdX[[HH[i]]] <- holdX[[HH[i]]]*0
	holdX[[HL[i]]] <- holdX[[HL[i]]]*0
	holdX[[LH[i]]] <- holdX[[LH[i]]]*0
	holdY[[HH[i]]] <- holdY[[HH[i]]]*0
        holdY[[HL[i]]] <- holdY[[HL[i]]]*0
        holdY[[LH[i]]] <- holdY[[LH[i]]]*0
	if(climate) {
	   holdClim[[HH[i]]] <- holdClim[[HH[i]]]*0
           holdClim[[HL[i]]] <- holdClim[[HL[i]]]*0
           holdClim[[LH[i]]] <- holdClim[[LH[i]]]*0
	}
      } # end of for 'i' loop.
    } # end of if do detail plots.
    par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)
    if(is.element("energy",which.plots)) {
	yl <- range(c(x$energy),na.rm=TRUE)
	plot(1:J, x$energy[,1], type="l", xlab="level", ylab="Energy", lwd=1.5, col="darkblue", ylim=yl, ...)
	lines(1:J, x$energy[,2], lty=2, col="darkblue", lwd=1.5, ylim=yl, ...)
	if(climate) {
	   lines(1:J, x$energy[,3], lty=3, col="green", lwd=1, ylim=yl, ...)
	   legend("topleft", legend=c(main1,main2,main3), col=c("darkblue", "darkblue", "green"), lty=1:3, lwd=c(1.5, 1.5, 1), bty="n")
	} else legend("topleft", legend=c(main1, main2), col="darkblue", lty=1:2, lwd=1.5, bty="n")
    } # end of plot energies stmts.
    if(any(is.element(c("all", "mse"),which.plots))) plot(1:J, x$mse, type="l", xlab="level", ylab="MSE", lwd=1.5, col="darkblue", ...)
    if(any(is.element(c("all", "rmse"),which.plots))) plot(1:J, x$rmse, type="l", xlab="level", ylab="RMSE", lwd=1.5, col="darkblue", ...)
    if(any(is.element(c("all", "acc"),which.plots)) & climate) plot(1:J, x$acc, type="l", xlab="level", ylab="ACC", lwd=1.5, col="darkblue", ...)
    if(!is.null(a$msg)) mtext(a$msg, line=0.05, outer=TRUE)
    invisible()
} # end of 'plot.waverify2d' function.

mowaverify2d <- function(X, ..., Clim=NULL, wavelet.type="haar", J=4) {
    UseMethod("mowaverify2d",X)
} # end of 'mowaverify2d' function.

mowaverify2d.SpatialVx <- function(X, ..., Clim=NULL, wavelet.type="haar", J=4,
				    useLL=FALSE, compute.shannon=FALSE, which.space="field",
				    time.point=1, model=1, verbose=FALSE) {

    out <- list()
    a <- attributes(X)
    attributes(out) <- a

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(X, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(X, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(X, model=model)
    else dat <- datagrabber(X)

    x <- dat$X
    y <- dat$Xhat
    ## End: Get the data sets

    out2 <- mowaverify2d.default(X=x, Clim=Clim, Y=y, wavelet.type=wavelet.type, J=J,
				    useLL=useLL, compute.shannon=compute.shannon, which.space=which.space, verbose=verbose)

    out$J <- out2$J
    out$X.wave <- out2$X.wave
    out$Y.wave <- out2$Y.wave
    out$energy <- out2$energy
    out$mse <- out2$mse
    out$rmse <- out2$rmse
    if(!is.null(out2$acc)) out$acc <- out2$acc
    out$wave.type <- out2$wave.type
    out$wave.method <- out2$wave.method

    attr(out, "time.point") <- time.point
    attr(out, "model") <- model
    class(out) <- "waverify2d"
    return(out)
} # end of 'mowaverify2d.SpatialVx' function.

mowaverify2d.default <- function(X, ..., Clim=NULL, Y, wavelet.type="haar", J=4,
				    useLL=FALSE, compute.shannon=FALSE, which.space="field", verbose=FALSE) {
   out <- list()
   attributes(out) <- atmp <- list(...)
   if(is.null(atmp$loc.byrow)) attr(out, "loc.byrow") <- FALSE

   out$J <- J
   climate <- !is.null(Clim)
   xdim <- dim(X)
   if(verbose) cat("\n", "Finding MO wavelet decompositions.\n")
   wv.X <- modwt.2d(X, wf=wavelet.type, J=J)
   wv.Y <- modwt.2d(Y, wf=wavelet.type, J=J)
   if( climate) wv.Clim <- modwt.2d( Clim, wf=wavelet.type, J=J)
   if(verbose) cat("\n", "MO wavelet decompositions found.\n")
   out$X.wave <- wv.X
   out$Y.wave <- wv.Y
   if( climate) out$Clim.wave <- wv.Clim
   lnames <- makeWaveNames(J=J)
   HH <- lnames$HH
   HL <- lnames$HL
   LH <- lnames$LH
   LL <- lnames$LL

   mse.i <- rmse.i <- numeric(J) + NA
   if( climate) {
	acc.i <- rmse.i
	energy.i <- matrix( NA, nrow=J, ncol=3)
	colnames(energy.i) <- c("Observed", "Forecast", "Climatology")
   } else {
	energy.i <- matrix( NA, nrow=J, ncol=2)
	colnames( energy.i) <- c("Observed", "Forecast")
   }

   if( compute.shannon) {
      if(verbose) cat("\n", "Cacluating Shannon entropy for each field.\n")
      Shannon.entropy <- matrix( NA, nrow=1, ncol=2)
      colnames( Shannon.entropy) <- c("Verification", "Forecast")
      tmp <- abs(unlist(wv.X))
      tmp <- tmp/sum(tmp,na.rm=TRUE)
      id0 <- tmp==0
      tmp[id0] <- 0
      tmp[!id0] <- tmp[!id0]*log(tmp[!id0])
      Shannon.entropy[,1] <- -sum(tmp, na.rm=TRUE)
      tmp <- abs(unlist(wv.Y))
      tmp <- tmp/sum(tmp,na.rm=TRUE)
      id0 <- tmp==0
      tmp[id0] <- 0
      tmp[!id0] <- tmp[!id0]*log(tmp[!id0])
      Shannon.entropy[,2] <- -sum(tmp, na.rm=TRUE)
      out$Shannon.entropy <- Shannon.entropy
      if(verbose) cat("\n", "Shannon entropy found.\n")
   }

   if(verbose) cat("n", "Calculating statistics for each detail field.\n")
   for(i in 1:J) {
	energy.i[i,1] <- sum( wv.X[[HH[i]]]^2, na.rm=TRUE) + sum( wv.X[[HL[i]]]^2, na.rm=TRUE) + sum( wv.X[[LH[i]]]^2, na.rm=TRUE)
	energy.i[i,2] <- sum( wv.Y[[HH[i]]]^2, na.rm=TRUE) + sum( wv.Y[[HL[i]]]^2, na.rm=TRUE) + sum( wv.Y[[LH[i]]]^2, na.rm=TRUE)
	if( climate) energy.i[i,3] <- sum( wv.Clim[[HH[i]]]^2, na.rm=TRUE) + sum( wv.Clim[[HL[i]]]^2, na.rm=TRUE) + sum( wv.Clim[[LH[i]]]^2, na.rm=TRUE)

	invX <- detailer(wv.X, level=i, which.space=which.space, trans.type="MODWT", useLL=useLL, lnames=lnames, J=J)
	invY <- detailer(wv.Y, level=i, which.space=which.space, trans.type="MODWT", useLL=useLL, lnames=lnames, J=J)
	if( climate) invClim <- detailer(wv.Clim, level=i, which.space=which.space, trans.type="MODWT", useLL=useLL, lnames=lnames, J=J)

	mse.i[i] <- mean( (c(invX) - c(invY))^2, na.rm=TRUE)
   	rmse.i[i] <- sqrt( mse.i[i])

	if( climate) {
           denom <- sqrt( sum( colSums((c( invX) - c( invClim))^2,na.rm=TRUE), na.rm=TRUE))*sqrt( sum( colSums((c(invY) - c(invClim))^2,na.rm=TRUE), na.rm=TRUE))
           numer <- sum( diag( t(invX - invClim) %*% (invY - invClim)), na.rm=TRUE)
           acc.i[i] <- numer/denom
	}

   } # end of for 'i' loop.
   out$energy <- energy.i
   out$mse <- mse.i
   out$rmse <- rmse.i
   if(climate) out$acc <- acc.i
   out$wave.type <- wavelet.type
   out$wave.method <- "MODWT"
   class(out) <- "waverify2d"
   return(out)
} # end of 'mowaverify2d.default' function.

waveIS <- function(object, J=NULL, wavelet.type="haar", time.point = 1, model = 1, 
    levels = NULL, max.n = NULL, smooth.fun = "hoods2dsmooth", 
    smooth.params = NULL, verbose=FALSE) {

    object <- hoods2dPrep(object=object, levels=levels, max.n=max.n,
		    smooth.fun=smooth.fun, smooth.params=smooth.params) 
   
    object.attr <- attributes(object)
     
   if (verbose) begin.time <- Sys.time()
   out <- list()
   attributes(out) <- object.attr

   thresholds <- object.attr$thresholds
   q <- dim(thresholds)[1]

   ## Begin: Get the data sets
   if(!missing(time.point) && !missing(model)) dat <- datagrabber(object, time.point=time.point, model=model)
   else if(!missing(time.point)) dat <- datagrabber(object, time.point=time.point)
   else if(!missing(model)) dat <- datagrabber(object, model=model)
   else dat <- datagrabber(object)

   X <- dat$X
   Y <- dat$Xhat
   ## End: Get the data sets

   xdim <- object.attr$xdim
   bigN <- prod(xdim)
   binmat <- matrix(0, xdim[1], xdim[2])
   Jtry <- log2(xdim)
   if(all(floor(Jtry) == ceiling(Jtry))) dyadic <- TRUE
   else dyadic <- FALSE
   if(is.null(J)) {
	if(dyadic) J <- min(Jtry,na.rm=TRUE)
	else J <- 4
   }
   if(dyadic) wmeth <- "DWT"
   else wmeth <- "MODWT"
   out$wave.method <- wmeth

   eX <- eY <- MSE <- SS <- matrix(NA, nrow=J, ncol=q)
   Bias <- mserand <- numeric(q)+NA

   if(verbose) cat("\n", "Looping through thresholds = \n")
   for(threshold in 1:q) {
	if(verbose) cat(thresholds[threshold], " ")
	Xbin <- Ybin <- binmat
	Ybin[Y>=thresholds[threshold,1]] <- 1
	Xbin[X>=thresholds[threshold,2]] <- 1
	s <- sum(colSums(Xbin,na.rm=TRUE),na.rm=TRUE)
        B <- sum(colSums(Ybin,na.rm=TRUE),na.rm=TRUE)/s
        s <- s/bigN
        Bias[threshold] <- B
        MSE.random <- B*s*(1-s) + s*(1-B*s)
	mserand[threshold] <- MSE.random
	if(dyadic) {
	   wv.X <- dwt.2d(Xbin,wf=wavelet.type,J=J)
	   wv.Y <- dwt.2d(Ybin,wf=wavelet.type,J=J)
	   wv.diff <- dwt.2d(Ybin - Xbin,wf=wavelet.type,J=J)
	} else {
	   wv.X <- modwt.2d(Xbin,wf=wavelet.type,J=J)
	   wv.Y <- modwt.2d(Ybin,wf=wavelet.type,J=J)
	   wv.diff <- modwt.2d(Ybin - Xbin,wf=wavelet.type,J=J)
	}
	if(threshold==1) lnames <- makeWaveNames(J=J)
	eX[,threshold] <- energizer(wv.X,lnames=lnames,J=J)
	eY[,threshold] <- energizer(wv.Y,lnames=lnames,J=J)
	for(j in 1:J) {
	   hold <- detailer(wv.diff,level=j,which.space="wavelet",trans.type=wmeth,lnames=lnames,J=J)
	   N <- length(hold[!is.na(hold)])
	   mse.ui <- sum(colSums(hold,na.rm=TRUE),na.rm=TRUE)/N
	   MSE[j,threshold] <- mse.ui
	   SS[j,threshold] <- 1 - (mse.ui*(J+1))/MSE.random
	} # end of for 'j' loop.
   } # end of for 'threshold' loop.
   if(verbose) cat("\n", "Finished computing DWTs for each threshold.\n")
   out$J <- J
   out$EnVx <- eX
   out$EnFcst <- eY
   out$MSE <- MSE
   out$Bias <- Bias
   out$SS <- SS
   out$MSE.random <- mserand
   class(out) <- "waveIS"
   return(out)
} # end of 'waveIS' function.

energizer <- function(x,lnames=NULL,J=NULL) {
   if(is.null(J)) J <- (length(x)-1)/3
   if(is.null(lnames)) lnames <- makeWaveNames(J=J)
   HH <- lnames$HH
   HL <- lnames$HL
   LH <- lnames$LH
   LL <- lnames$LL
   out <- numeric(J)+NA
   for(i in 1:J) {
	good <- !(is.na(x[[HH[i]]]) & is.na(x[[HL[i]]]) & is.na(x[[LH[i]]]))
	N <- length(x[[HH[i]]][good])
	out[i] <- sum(colSums(x[[HH[i]]]^2 + x[[HL[i]]]^2 + x[[LH[i]]]^2,na.rm=TRUE), na.rm=TRUE)/N
   } # end of for 'i' loop.
   return(out)
} # end of 'energizer' function.

summary.waveIS <- function(object,...) {
   x <- object
   a <- attributes(object)
   args <- list(...)
   if(!is.null(args$silent)) silent <- args$silent
   else silent <- FALSE
   J <- x$J
   MSE <- x$MSE
   colnames(MSE) <- a$qs
   MSEu <- apply(MSE,2,sum,na.rm=TRUE)
   MSEperc <- t(t(MSE)/MSEu)*100
   SSu <- matrix(1 - MSEu/x$MSE.random, nrow=1)
   SS <- x$SS
   colnames(SS) <- a$qs
   EnVx <- x$EnVx
   colnames(EnVx) <- a$qs
   EnFcst <- x$EnFcst
   colnames(EnFcst) <- a$qs
   EnVx.u <- apply(EnVx,2,sum,na.rm=TRUE)
   EnFcst.u  <- apply(EnFcst,2,sum,na.rm=TRUE)
   EnVx.perc <- t(t(EnVx)/EnVx.u)*100
   EnFcst.perc <- t(t(EnFcst)/EnFcst.u)*100
   EnRelDiff <- (EnFcst - EnVx)/(EnFcst + EnVx)
   Bias <- matrix(x$Bias, nrow=1)
   colnames(Bias) <- a$qs

   xdim <- a$xdim
   dyadic <- all(floor(log2(xdim))==ceiling(log2(xdim)))
   if(dyadic) level.names <- paste("2^",J-(1:J),sep="")
   else level.names <- as.character(1:J)

   rownames(MSE) <- rownames(EnVx) <- rownames(SS) <- rownames(MSEperc) <- 
	rownames(EnFcst) <- rownames(EnVx.perc) <- rownames(EnFcst.perc) <- 
	rownames(EnRelDiff) <- level.names

   if(!silent) {
	cat("\n", x$Fcst.name, " model field compared with verification field ", x$Vx.name, "\n")
	cat("MSE by threshold (columns) and level (rows).\n")
	print(round(MSE, digits=4))
	cat("\n", "MSE% by threshold (columns) and level (rows).\n")
	print(round(MSEperc, digits=0))
	cat("\n", "Field skill score by threshold.\n")
	print(round(SSu, digits=4))
	cat("\n", "IS Skill Score\n")
	print(round(SS, digits=4))
	cat("\n", "Verification energy\n")
	print(round(EnVx, digits=4))
	cat("\n", "% Verification energy\n")
	print(round(EnVx.perc,digits=0))
	cat("\n", "Forecast energy\n")
	print(round(EnFcst, digits=4))
	cat("\n", "% Forecast energy\n")
	print(round(EnFcst.perc, digits=0))
	cat("\n", "Energy Relative Difference\n")
	print(round(EnRelDiff, digits=4))
	cat("\n", "Frequency Bias\n")
	print(round(Bias, digits=4))
   } # end of if '!silent' stmts.
   invisible(c(x,list(MSEu=MSEu,MSEperc=MSEperc,SSu=SSu,EnVx.u=EnVx.u,EnFcst.u=EnFcst.u,EnVx.perc=EnVx.perc,EnFcst.perc=EnFcst.perc,EnRelDiff=EnRelDiff)))
} # end of 'summary.waveIS' function.

plot.waveIS <- function(x, main1="X", main2="Y", which.plots=c("all", "mse", "ss", "energy"), level.label=NULL, ...) {
   a <- attributes(x)
   x <- summary(x,silent=TRUE)
   J <- x$J
   xdim <- a$xdim
   thresholds <- a$qs
   q <- length(thresholds)

#   if(!is.null(a$levels)) {
#	levels <- a$levels
#	if(is.null(level.label)) level.label <- "scale/level"
#   } else {
	dyadic <- all(floor(log2(xdim))==ceiling(log2(xdim)))
	if(dyadic) levels <- paste("2^", J-(1:J), sep="")
	else {
	   levels <- 1:J
	   if(is.null(level.label)) level.label <- "Level"
	}
#   }

   if(any(is.element(c("all","mse"),which.plots))) {
	par(mfrow=c(2,1), mar=c(5.1,4.1,4.1,5.1))
	image(x$MSE, col=c("grey",tim.colors(64)), main=paste("MSE (", main1, " vs ", main2, ")", sep=""), xlab=level.label,
		ylab=paste("threshold (", a$units, ")", sep=""), axes=FALSE)
	axis(1, at=seq(0,1,,J), labels=levels)
	axis(2, at=seq(0,1,,q), labels=thresholds)
	image.plot(x$MSE, col=c("grey",tim.colors(64)),legend.only=TRUE)
	MSEperc <- x$MSEperc
	image(MSEperc, col=c("grey",tim.colors(64)), zlim=c(0,100), main="MSE %", xlab=level.label,
				ylab=paste("threshold (", a$units, ")", sep=""), axes=FALSE)
        axis(1, at=seq(0,1,,J), labels=levels)
        axis(2, at=seq(0,1,,q), labels=thresholds)
        image.plot(MSEperc, col=c("grey",tim.colors(64)), zlim=c(0,100),legend.only=TRUE)
   } # end of if do MSE stmts.
   if(any(is.element(c("all","ss"),which.plots))) {
	par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,5.1))
	image(x$SS, col=c("grey",tim.colors(64)), main=paste("IS Skill Score (", main1, " vs ", main2, ")", sep=""), xlab=level.label,
				ylab=paste("threshold (", a$units, ")", sep=""), axes=FALSE)
        axis(1, at=seq(0,1,,J), labels=levels)
        axis(2, at=seq(0,1,,q), labels=thresholds)
        image.plot(x$SS, col=c("grey",tim.colors(64)),legend.only=TRUE)
   } # end of if do SS stmts.
   if(any(is.element(c("all","energy"),which.plots))) {
	zl <- range(c(c(x$EnVx),c(x$EnFcst)),finite=TRUE)
	par(mfrow=c(3,2), mar=c(5.1,4.1,4.1,5.1))
	image(x$EnVx, col=c("grey",tim.colors(64)), main=paste(main1, " Energy", sep=""), xlab=level.label,
				ylab=paste("threshold (", a$units, ")", sep=""), axes=FALSE, zlim=zl)
        axis(1, at=seq(0,1,,J), labels=levels)
        axis(2, at=seq(0,1,,q), labels=thresholds)
        image.plot(x$EnVx, col=c("grey",tim.colors(64)),legend.only=TRUE, zlim=zl)
        image(x$EnFcst, col=c("grey",tim.colors(64)), zlim=zl, main=paste(main2, " Energy", sep=""), xlab=level.label,
                                ylab=paste("threshold (", a$units, ")", sep=""), axes=FALSE)
        axis(1, at=seq(0,1,,J), labels=levels)
        axis(2, at=seq(0,1,,q), labels=thresholds)
        image.plot(x$EnFcst, col=c("grey",tim.colors(64)), zlim=zl,legend.only=TRUE)

	image(x$EnVx.perc, col=c("grey",tim.colors(64)), main=paste(main1, " Energy %", sep=""), xlab=level.label,
                                ylab=paste("threshold (", a$units, ")", sep=""), axes=FALSE, zlim=c(0,100))
        axis(1, at=seq(0,1,,J), labels=levels)
        axis(2, at=seq(0,1,,q), labels=thresholds)
        image.plot(x$EnVx.perc, col=c("grey",tim.colors(64)),legend.only=TRUE, zlim=c(0,100))
        image(x$EnFcst.perc, col=c("grey",tim.colors(64)), zlim=c(0,100), main=paste(main2, " Energy %", sep=""), xlab=level.label,
                                ylab=paste("threshold (", a$units, ")", sep=""), axes=FALSE)
        axis(1, at=seq(0,1,,J), labels=levels)
        axis(2, at=seq(0,1,,q), labels=thresholds)
        image.plot(x$EnFcst.perc, col=c("grey",tim.colors(64)), zlim=c(0,100),legend.only=TRUE)

	image(x$EnRelDiff, col=c("grey",tim.colors(64)), main=paste("Energy Relative Difference (", main1, " vs ", main2, ")", sep=""), xlab=level.label,
                                ylab=paste("threshold (", a$units, ")", sep=""), axes=FALSE, zlim=c(-1,1))
        axis(1, at=seq(0,1,,J), labels=levels)
        axis(2, at=seq(0,1,,q), labels=thresholds)
        image.plot(x$EnRelDiff, col=c("grey",tim.colors(64)), zlim=c(-1,1),legend.only=TRUE)
   } # end of if plot energy results.
   invisible()
} # end of 'plot.waveIS' function.

detailer <- function(x, level, which.space="field", trans.type="DWT", useLL=FALSE, lnames=NULL, J=NULL) {
   if(is.null(lnames)) {
	if(is.null(J)) stop("detailer: must supply J if lnames is NULL.")
	lnames <- makeWaveNames(J=J)
   }
   HH <- lnames$HH
   HL <- lnames$HL
   LH <- lnames$LH
   LL <- lnames$LL

   if(which.space=="field") {
	hold <- makeEmptyWave(x)
	if(useLL) hold[[LL]] <- x[[LL]]
	hold[[HH[level]]] <- x[[HH[level]]]
	hold[[HL[level]]] <- x[[HL[level]]]
	hold[[LH[level]]] <- x[[LH[level]]]
	ni <- length(x[[HH[level]]])*3
	lam <- sqrt(sum(c(hold[[HH[level]]]^2)+c(hold[[HL[level]]]^2)+c(hold[[LH[level]]]^2),na.rm=TRUE)/ni)
	hold[[HH[level]]] <- sign(hold[[HH[level]]])*pmax(c(abs(hold[[HH[level]]])-lam),0)
	hold[[HL[level]]] <- sign(hold[[HL[level]]])*pmax(c(abs(hold[[HL[level]]])-lam),0)
	hold[[LH[level]]] <- sign(hold[[LH[level]]])*pmax(c(abs(hold[[LH[level]]])-lam),0)
	if(trans.type=="DWT") out <- idwt.2d(hold)
	else out <- imodwt.2d(hold)
   } else out <- (x[[HH[level]]]/(2^level))^2 + (x[[HL[level]]]/(2^level))^2 + (x[[LH[level]]]/(2^level))^2
   return(out)
} # end of 'detailer' function.

wavePurifyVx <- function(x,y,object=NULL, climate=NULL,
    which.stats=c("bias", "ts", "ets", "pod", "far", "f", "hk", "mse"),
    thresholds=NULL, return.fields=FALSE, time.point=1, model=1, verbose=FALSE, ...) {

   if(verbose) begin.time <- Sys.time()
   isnt.xy <- missing(x) | missing(y)
   if(isnt.xy & !is.null(object)) {
	out <- list()
	out.attr <- attributes(object)
	attributes(out) <- out.attr

	## Begin: Get the data sets
        if(!missing(time.point) && !missing(model)) dat <- datagrabber(object, time.point=time.point, model=model)
        else if(!missing(time.point)) dat <- datagrabber(object, time.point=time.point)
        else if(!missing(model)) dat <- datagrabber(object, model=model)
        else dat <- datagrabber(object)

        x <- dat$X
        y <- dat$Xhat
        ## End: Get the data sets

	xdim <- out.attr$xdim
	thresholds <- out.attr$thresholds

   } else if(!isnt.xy) {
	out <- list()
	attributes(out) <- atmp <- list(...)
	if(is.null(atmp$loc.byrow)) attr(out, "loc.byrow") <- FALSE
	xdim <- dim(x)
	if(is.null(thresholds) & any(is.element(c("bias", "ts", "ets", "pod", "far", "f", "hk"), which.stats))) {
	   thresholds <- cbind( quantile(c(x), probs=c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9, 0.95)),
				quantile(c(y), probs=c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9, 0.95)))
	   out$qs <- c("0", "0.1", "0.25", "0.33", "0.5", "0.66", "0.75", "0.9", "0.95")
	} else if(!is.matrix(thresholds)) thresholds <- cbind(thresholds,thresholds)
	colnames(thresholds) <- c("Forecast","Verification")
	out$thresholds <- thresholds
	out$Vx.name <- as.character(substitute(x))
	out$Fcst.name <- as.character(substitute(y))
   } else stop("wavePurifyVx: Must supply arguments x and y, or object.")

   args <- list(...)
   if(isnt.J <- is.null(args$J)) {
	J <- floor(log2(min(xdim)))
	args$J <- J
   } else {
	J <- args$J
   }
   out$args <- args
   if(all(floor(log2(xdim))==ceiling(log2(xdim)))) dyadic <- TRUE
   else dyadic <- FALSE
   if(verbose) cat("\n", "Denoising the fields.\n")
   if(dyadic) {
	if(isnt.J) {
	   X <- denoise.dwt.2d(x=x, J=J, ...)
	   Y <- denoise.dwt.2d(x=y, J=J, ...)
	   if(!is.null(climate)) Climate <- denoise.dwt.2d(x=climate, J=J, ...)
	} else {
	   X <- denoise.dwt.2d(x=x, ...)
	   Y <- denoise.dwt.2d(x=y, ...)
	   if(!is.null(climate)) Climate <- denoise.dwt.2d(x=climate, ...)
	} # end of if else 'isnt.J' stmts.
   } else {
	if(isnt.J) {
           X <- denoise.modwt.2d(x=x, J=J, ...)
           Y <- denoise.modwt.2d(x=y, J=J, ...)
	   if(!is.null(climate)) Climate <- denoise.modwt.2d(x=climate, J=J, ...)
        } else {
           X <- denoise.modwt.2d(x=x, ...)
           Y <- denoise.modwt.2d(x=y, ...)
	   if(!is.null(climate)) Climate <- denoise.modwt.2d(x=climate, ...)
        } # end of if 'is.J' stmts.
   } # end of if else 'dyadic' stmts.
   if(return.fields) {
	out$X <- x
	out$Xhat <- y
	out$X2 <- X
	out$Y2 <- Y
	if(!is.null(climate)) {
	   out$Climate <- climate
	   out$Climate2 <- Climate
	}
   }
   if(!is.null(thresholds)) q <- dim(thresholds)[1]
   else q <- 1
   if(is.element("bias",which.stats)) out$bias <- numeric(q)+NA 
   if(is.element("ts",which.stats)) out$ts <- numeric(q)+NA 
   if(is.element("ets",which.stats)) out$ets <- numeric(q)+NA 
   if(is.element("pod",which.stats)) out$pod <- numeric(q)+NA 
   if(is.element("far",which.stats)) out$far <- numeric(q)+NA 
   if(is.element("f",which.stats)) out$f <- numeric(q)+NA 
   if(is.element("hk",which.stats)) out$hk <- numeric(q)+NA 
   if(is.element("mse",which.stats)) out$mse <- numeric(q)+NA
   if(!is.null(climate)) out$acc <- numeric(q)+NA
   binmat <- matrix(0, xdim[1], xdim[2])
   if(verbose) {
	cat("\n", "Looping through thresholds = \n") 
	print(thresholds)
	cat("\n")
   }
   for( threshold in 1:q) {
      if(verbose) cat(threshold, " ")
      if(is.element("mse",which.stats)) {
	if(!is.null(thresholds)) {
	   X2 <- X
	   Y2 <- Y
	   X2[X<thresholds[threshold,2]] <- 0
	   Y2[Y<thresholds[threshold,1]] <- 0
	   out$mse[threshold] <- vxstats(Y2,X2,which.stats="mse")$mse
	} else out$mse[threshold] <- vxstats(Y,X,which.stats="mse")$mse
      } # end of if do mse stmts.
      if(any(is.element(c("bias", "ts", "ets", "pod", "far", "f", "hk"), which.stats))) {
	Xbin <- Ybin <- binmat
	Xbin[X>=thresholds[threshold,2]] <- 1
	Ybin[Y>=thresholds[threshold,1]] <- 1
	if(threshold==1) {
	   dostats <- which.stats
	   if(is.element("mse",dostats)) dostats <- dostats[dostats != "mse"]
	}
	tmp <- vxstats(Ybin,Xbin,which.stats=dostats)
	if(is.element("bias",which.stats)) out$bias[threshold] <- tmp$bias 
	if(is.element("ts",which.stats)) out$ts[threshold] <- tmp$ts 
	if(is.element("ets",which.stats)) out$ets[threshold] <- tmp$ets 
	if(is.element("pod",which.stats)) out$pod[threshold] <- tmp$pod 
	if(is.element("far",which.stats)) out$far[threshold] <- tmp$far 
	if(is.element("f",which.stats)) out$f[threshold] <- tmp$f 
	if(is.element("hk",which.stats)) out$hk[threshold] <- tmp$hk
      } # end of if any vxstats besides mse stmts.
      if(!is.null(climate)) {
	if(!is.null(thresholds)) {
           X2 <- X
           Y2 <- Y
	   Clim <- Climate
           X2[X<thresholds[threshold,2]] <- 0
           Y2[Y<thresholds[threshold,1]] <- 0
	   Clim[climate<thresholds[threshold,1]] <- 0
	   denom <- sqrt(sum(colSums((Y2-Clim)^2,na.rm=TRUE),na.rm=TRUE))*sqrt(sum(colSums((X2-Clim)^2,na.rm=TRUE),na.rm=TRUE))
           numer <- sum(diag(t(Y2-Clim) %*% (X2-Clim)),na.rm=TRUE)
        } else {
	   denom <- sqrt(sum(colSums((Y-Climate)^2,na.rm=TRUE),na.rm=TRUE))*sqrt(sum(colSums((X-Climate)^2,na.rm=TRUE),na.rm=TRUE))
	   numer <- sum(diag(t(Y-Climate) %*% (X-Climate)),na.rm=TRUE)
	}
	out$acc[threshold] <- numer/denom
      } # end of if do acc stmts.
   } # end of for 'threshold' loop.

   attr(out, "time.point") <- time.point
   attr(out, "model") <- model   
   class(out) <- "wavePurifyVx"
   if(verbose) print(Sys.time() - begin.time)
   return(out)
} # end of 'wavePurifyVx' function.

summary.wavePurifyVx <- function(object, ...) {

   object.attr <- attributes(object)

   if(!is.null(object.attr$thresholds)) {
	u <- Thresholds <- object.attr$thresholds
   	q <- dim(u)[1]
   	cat("\n", "Thresholds applied are:\n")
   	print(Thresholds)
	   if(is.null(object.attr$qs)) {
		if(all(u[,1]==u[,2])) ulab <- as.character(u[,1])
		else ulab <- as.character(1:q)
	   } else ulab <- object.attr$qs
	   if(!is.null(object$bias)) {
		bias <- matrix(object$bias, nrow=1)
		colnames(bias) <- ulab
		cat("\n", "Frequency Bias: \n")
		print(bias)
	   }
	   if(!is.null(object$ts)) {
		ts <- matrix(object$ts, nrow=1)
		colnames(ts) <- ulab
		cat("\n", "Threat Score: \n")
		print(ts)
	   }
	   if(!is.null(object$ets))  {
		ets <- matrix(object$ets, nrow=1)
		colnames(ets) <- ulab
		cat("\n", "Gilbert Skill Score: \n")
		print(ets)
	   }
	   if(!is.null(object$pod))  {
		pod <- matrix(object$pod, nrow=1)
		colnames(pod) <- ulab
		cat("\n", "Probability of Detecting an Event (POD): \n")
		print(pod)
	   }
	   if(!is.null(object$far)) {
		far <- matrix(object$far, nrow=1)
		colnames(far) <- ulab
		cat("\n", "False Alarm Ratio: \n")
		print(far)
	   }
	   if(!is.null(object$f)) {
		f <- matrix(object$f, nrow=1)
		colnames(f) <- ulab
		cat("\n", "False Alarm Rate: \n")
		print(f)
	   }
	   if(!is.null(object$hk)) {
		hk <- matrix(object$hk, nrow=1)
		colnames(hk) <- ulab
		cat("\n", "Hanssen-Kuipers Score: \n")
		print(hk)
	   }
	   if(!is.null(object$mse)) {
		mse <- matrix(object$mse, nrow=1)
		colnames(mse) <- ulab
		cat("\n", "MSE: \n")
		print(mse)
	   }
	   if(!is.null(object$acc)) {
		acc <- matrix(object$acc, nrow=1)
		colnames(acc) <- ulab
		cat("\n", "ACC: \n")
		print(acc)
	   }
   } else {
	if(!is.null(object$mse)) {
                mse <- matrix(object$mse, nrow=1)
		cat("\n", "MSE: \n")
                print(mse)
           }
           if(!is.null(object$acc)) {
                acc <- matrix(object$acc, nrow=1)
		cat("\n", "ACC: \n")
                print(acc)
           } 
   } # end of whether or not thresholds applied stmts.
   invisible()
} # end of 'summary.wavePurifyVx' function.

plot.wavePurifyVx <- function(x, ..., set.pw=FALSE, col, zlim, horizontal=TRUE, type=c("stats", "fields")) {

   type <- tolower(type)
   type <- match.arg(type)

   if(missing(col)) col <- c("gray", tim.colors(64))

   a <- attributes(x)
   loc.byrow <- a$loc.byrow
   nm <- names(x)

   par(oma=c(0,0,2,0))

   if(!is.logical(set.pw) && length(set.pw) == 4 && is.numeric(set.pw)) {
	par(mfrow=set.pw)
	set.pw <- FALSE
   } else if(!is.logical(set.pw)) stop("plot.wavePurifyVx: invalid set.pw argument.")

   denX <- !is.null(x$X2)

   if(set.pw) {
        if(type=="fields" && denX) numplots <- 4
        else numplots <- 1
        if(type=="fields" && denX && !is.null(x$Climate)) numplots <- numplots+2
        if(type=="stats") {
            if(is.element("mse",nm) && any(is.element(c("bias","ts","ets","pod","far","hk"),nm))) numplots <- numplots+1
            if(is.element("acc",nm) && (is.element("mse",nm) | any(is.element(c("bias","ts","ets","pod","far","hk"),nm)))) numplots <- numplots+1
	}
        if(numplots > 1) {
	    if(is.null(x$Climate) && set.pw) par(mfrow=c(ceiling(numplots/2),2))
            else par(mfrow=c(ceiling(numplots/3),3))
        }
   }

   if(denX && type=="fields") {
        domap <- a$map
	if(is.null(domap)) domap <- FALSE

	loc <- a$loc

	if(domap) {
	    locr <- apply(loc,2,range,finite=TRUE)
	    ax <- list(x=pretty(round(loc[,1], digits=2)), y=pretty(round(loc[,2], digits=2)))
	}

	X <- x$X
	Xhat <- x$Xhat
	X2 <- x$X2
	Y2 <- x$Y2

	xd <- a$xdim
	if(is.null(xd)) xd <- dim(X)
	if(is.null(loc)) loc <- cbind(rep(1:xd[1],xd[2]), rep(1:xd[2],each=xd[1]))

	ptitle <- a$data.name
	if(is.null(ptitle)) ptitle <- c("","")
	else {
	    if(length(ptitle) == a$nforecast + 2) ptitle <- ptitle[-1]

	    if(a$nforecast > 1) {
	        if(!is.numeric(a$model)) {
		    modnum <- (1:length(ptitle))[ptitle==a$model] - 1
	        } else modnum <- a$model
	            ptitle <- c(ptitle[1], ptitle[modnum + 1])
	    }

	    tiid <- a$time
	    if(length(tiid) > 1) ptitle <- paste(ptitle, " (t = ", tiid[a$time.point], ")", sep="")
	} # end of if 'ptitle' is NULL stmts.

	if(!is.null(x$Climate)) {
		Clim <- x$Climate
		Clim2 <- x$Climate2
		if(missing(zlim)) zl <- range( c(c(X), c(Xhat), c(X2), c(Y2), c(Clim), c(Clim2)), finite=TRUE)
	} else if(missing(zlim)) zl <- range( c(c(X), c(Xhat), c(X2), c(Y2)), finite=TRUE)

	proj <- a$projection
	if(is.null(proj)) proj <- FALSE

	if(proj) {
	    if(domap) {
		map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
		poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow), matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
		    X, main=ptitle[1], col=col, zlim=zl, add=TRUE)
		map(add=TRUE, lwd=1.5)
		map(add=TRUE, database="state")

		map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
                poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow), matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                    Xhat, main=ptitle[2], col=col, zlim=zl, add=TRUE)
                map(add=TRUE, lwd=1.5)
                map(add=TRUE, database="state")

		if(!is.null(x$Climate)) {
                    map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
                    poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow), matrix(loc[,2], xd[1], xd[2],
                        byrow=loc.byrow), Clim, main="Climate", col=col, zlim=zl, add=TRUE)
                    map(add=TRUE, lwd=1.5)
                    map(add=TRUE, database="state")
                }

		map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
                poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow), matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                    X2, main=paste(ptitle[1], " de-noised", sep=""), col=col, zlim=zl, add=TRUE)
                map(add=TRUE, lwd=1.5)
                map(add=TRUE, database="state")

		map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
                poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow), matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                    Y2, main=paste(ptitle[2], " de-noised", sep=""), col=col, zlim=zl, add=TRUE)
                map(add=TRUE, lwd=1.5)
                map(add=TRUE, database="state")

		if(!is.null(x$Climate)) {
                    map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
                    poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow), matrix(loc[,2], xd[1], xd[2],
                        byrow=loc.byrow), Clim2, main="Climate de-noised", col=col, zlim=zl, add=TRUE)
                    map(add=TRUE, lwd=1.5)
                    map(add=TRUE, database="state")
                }

	    } else {
		poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow), matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                    X, main=ptitle[1], col=col, zlim=zl)
		poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow), matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                    Xhat, main=ptitle[2], col=col, zlim=zl)
		if(!is.null(x$Climate)) {
		    poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow),
				matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow), Clim, main="Climate", col=col, zlim=zl)
		}
		poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow), matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                    X2, main=paste(ptitle[1], " de-noised", sep=""), col=col, zlim=zl)
		poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow), matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                    Y2, main=paste(ptitle[2], " de-noised", sep=""), col=col, zlim=zl)
		if(!is.null(x$Climate)) {
                    poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow),
                                matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow), Clim2, main="Climate de-noised",
				col=col, zlim=zl)
                }
	    }
	} else {
	    if(domap) {
                map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
                image(as.image(X, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main=ptitle[1], col=col, zlim=zl, add=TRUE)
                map(add=TRUE, lwd=1.5)
                map(add=TRUE, database="state")

                map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
                image(as.image(Xhat, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main=ptitle[2], col=col, zlim=zl, add=TRUE)
                map(add=TRUE, lwd=1.5)
                map(add=TRUE, database="state")

                if(!is.null(x$Climate)) {
                    map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
                    image(as.image(Clim, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main="Climate", col=col, zlim=zl, add=TRUE)
                    map(add=TRUE, lwd=1.5)
                    map(add=TRUE, database="state")
                }

                map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
                image(as.image(X2, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main=paste(ptitle[1], " de-noised", sep=""),
		    col=col, zlim=zl, add=TRUE)
                map(add=TRUE, lwd=1.5)
                map(add=TRUE, database="state")

                map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
                image(as.image(Y2, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main=paste(ptitle[2], " de-noised", sep=""),
		    col=col, zlim=zl, add=TRUE)
                map(add=TRUE, lwd=1.5)
                map(add=TRUE, database="state")

                if(!is.null(x$Climate)) {
                    map(xlim=locr[,1], ylim=locr[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)
                    image(as.image(Clim2, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main="Climate de-noised", col=col, zlim=zl, add=TRUE)
                    map(add=TRUE, lwd=1.5)
                    map(add=TRUE, database="state")
                }

            } else {
                image(as.image(X, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main=ptitle[1], col=col, zlim=zl)
                image(as.image(Xhat, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main=ptitle[2], col=col, zlim=zl)
                if(!is.null(x$Climate)) {
                    image(as.image(Clim, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main="Climate", col=col, zlim=zl)
                }
                image(as.image(X2, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main=paste(ptitle[1], " de-noised", sep=""),
		    col=col, zlim=zl)
                image(as.image(Y2, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main=paste(ptitle[2], " de-noised", sep=""),
		    col=col, zlim=zl)
                if(!is.null(x$Climate)) {
                    image(as.image(Clim2, nx=xd[1], ny=xd[2], x=loc, na.rm=TRUE), main="Climate de-noised",
                                col=col, zlim=zl)
                }
            }
	} # end of if else 'projection' stmts.
        image.plot(X, zlim=zl, col=col, horizontal=horizontal, legend.only=TRUE)	
	return(invisible())
   } else if(!denX && type == "fields") stop("plot.wavePurifyVx: return.fields must be TRUE in call to wavePurifyVx to use type fields.")

    if(!is.null(a$thresholds)) {
      q <- dim(a$thresholds)[1]
      if(q  >= 2) {
	 if(is.null(a$units)) unt <- NULL
         else unt <- paste("(", a$units, ")", sep="")
         if(any(is.element(c("bias","ts","ets","pod","far","hk"),nm))) {
		if(is.element("hk",nm)) yl <- c(-1,1)
		else if(is.element("ets",nm)) yl <- c(-1/3,1)
		else yl <- c(0,1)
		plot(1:q, runif(q,-1,1), type="n", xaxt="n", ylim=yl, xlab=paste("Threshold ", unt, sep=""), ylab="")
		legnames <- c("TS", "GSS", "POD", "FAR", "F", "HK", "Bias")[is.element(c("ts", "ets", "pod", "far", "f", "hk", "bias"),nm)]
		nl <- length(legnames)
		if(nl==7) cl <- c(1:6,8)
		else cl <- 1:nl
		legend("topright", legend=legnames, lty=cl, col=cl, bty="n")
		i <- 1
		if(is.element("ts", nm)) {
		   lines(1:q, x$ts, col=i, lty=i, lwd=1.5)
		   i <- i+1
		}
		if(is.element("ets", nm)) {
		   lines(1:q, x$ets, col=i, lty=i, lwd=1.5)
                   i <- i+1
		}
		if(is.element("pod", nm)) {
                   lines(1:q, x$pod, col=i, lty=i, lwd=1.5)
                   i <- i+1
                }
		if(is.element("far", nm)) {
                   lines(1:q, x$far, col=i, lty=i, lwd=1.5)
                   i <- i+1
                }
		if(is.element("f", nm)) {
                   lines(1:q, x$f, col=i, lty=i, lwd=1.5)
                   i <- i+1
                }
		if(is.element("hk", nm)) {
                   lines(1:q, x$hk, col=i, lty=i, lwd=1.5)
                   i <- i+1
                }
		if(is.element("bias", nm)) {
		   if(i==7) i <- 8
		   brg <- range(x$bias, finite=TRUE)
		   if(length(unique(brg))==1) brg <- c(brg[1]-0.01, brg[2]+0.01)
		   par(usr=c(1,q,brg))
                   lines(1:q, x$bias, col=i, lty=i, lwd=1.5)
		   bax <- pretty(sort(unique(x$bias)))
		   axis(4, at=bax, labels=bax, col=i)
                }
      	   } # end of if 'any c("bias", "ts", "ets", "pod", "far", "f", "hk")' to make in plot stmts.
	   if(is.element("mse",nm)) plot(1:q, x$mse, type="l", col="darkblue", xlab=paste("Threshold ", unt, sep=""), ylab="MSE")
	   if(is.element("acc",nm)) plot(1:q, x$acc, type="l", col="darkblue", xlab=paste("Threshold ", unt, sep=""), ylab="ACC")
	} # end of if more than one threshold stmt.
   } # end of if 'thresholds' stmts.
   mtext(a$msg, line=0.05, outer=TRUE)
   invisible()
} # end of 'plot.wavePurifyVx' function.
