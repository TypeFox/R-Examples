initGMM <- function(x) {
   x <- disjointer(x)
   a <- lapply(x, FeatureProps, which.props=c("centroid", "area"))
   cent <- matrix( unlist(lapply(a, function(x) return(x$centroid))), ncol=2, byrow=TRUE)
   a <- unlist(lapply(a, function(x) return(x$area)))
   o <- order(a,decreasing=TRUE)
   ifun <- function(x) {
	y <- FeatureAxis(x)
	if(is.null(y$aspect.ratio)) return(c(1, 1))
	if(y$aspect.ratio>1) return(c(y$lengths$MinorAxis, y$lengths$MajorAxis))
	else return(c(y$lengths$MajorAxis, y$lengths$MinorAxis))
   } # end of internal 'ifun' function.
   s <- matrix( unlist(lapply(x, ifun)), ncol=2, byrow=TRUE)
   return(data.frame(ind=o, area=a, Xcen=cent[,1], Ycen=cent[,2], sX=s[,1], sY=s[,2]))
} # end of 'initGMM function.

gmm2d <- function(x, ...) {
    UseMethod("gmm2d", x)
} # end of 'gmm2d' function.

gmm2d.SpatialVx <- function(x, ..., time.point=1, model=1, K=3, gamma=1, threshold=NULL, initFUN="initGMM",
    verbose=FALSE) {

    a <- attributes(x)

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(x, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(x, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(x, model=model)
    else dat <- datagrabber(x)
   
    X <- dat$X
    Xhat <- dat$Xhat
    ## End: Get the data sets

    out <- gmm2d.default(x=X, xhat=Xhat, K=K, gamma=gamma, threshold=threshold, initFUN=initFUN, verbose=verbose)

    if(length(a$data.name) == a$nforecast + 2) {
        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[1:2]
    } else {
        dn <- a$data.name[-1]
        vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    attr(out, "data.name") <- c(vxname, dn[model.num])
    attr(out, "projection") <- a$projection
    attr(out, "map") <- a$map
    attr(out, "loc") <- a$loc
    attr(out, "loc.byrow") <- a$loc.byrow
    attr(out, "msg") <- a$msg

    return(out)
} # end of 'gmm2d.SpatialVx' function.

gmm2d.default <- function(x, ..., xhat, K=3, gamma=1, threshold=NULL, initFUN="initGMM", verbose=FALSE) {
    X <- x
   xdim <- dim(X)
   tmpX <- X
   if(!is.null(threshold)) tmpX[X<threshold] <- 0
   tmpY <- xhat
   if(!is.null(threshold)) tmpY[xhat<threshold] <- 0
   xdim <- dim(X)
   Ax <- sum(colSums(tmpX,na.rm=TRUE),na.rm=TRUE)
   Ay <- sum(colSums(tmpY,na.rm=TRUE),na.rm=TRUE)
   if(verbose) cat("\n", "Finding initial values for verification/validation field.\n")
   # x <- initGMM(tmpX)
   x <- do.call(initFUN, list(tmpX,...))
   if(verbose) cat("\n", "Finding initial values for forecast/model field.\n")
   # y <- initGMM(tmpY)
   y <- do.call(initFUN, list(tmpY,...))
   if(verbose) cat("\n", "Finding initial means for the first ", K, " largest identified objects.\n")
   meanX <- c(unlist(x[x$ind[1:K],3]), unlist(x[x$ind[1:K],4]))
   meanY <- c(unlist(y[y$ind[1:K],3]), unlist(y[y$ind[1:K],4]))
   sdX   <- rbind(c(unlist(x[x$ind[1:K],5])), rep(0,K), rep(0,K), c(unlist(x[x$ind[1:K],6])))
   sdY   <- rbind(c(unlist(y[y$ind[1:K],5])), rep(0,K), rep(0,K), c(unlist(y[y$ind[1:K],6])))
   # Find the value of m (Eq (11) from Lakshmanan and Kain 2010) for each field and non-zero grid point.
   if(verbose) cat("\n", "Setting up locations based on intensities.\n")
   sX <- sY <- cbind(rep(1:xdim[1],xdim[2]),rep(1:xdim[2],each=xdim[1]))
   tmpXc <- c(tmpX)
   tmpYc <- c(tmpY)
   sX <- sX[tmpXc>0,]
   sY <- sY[tmpYc>0,]
   df.sX <- hist(tmpXc[tmpXc>0],plot=FALSE)
   df.sY <- hist(tmpYc[tmpYc>0],plot=FALSE)
   nX <- length(df.sX$breaks)
   nY <- length(df.sY$breaks)
   bX <- df.sX$breaks[2:nX]
   bY <- df.sY$breaks[2:nY]
   cdfX <- cumsum(df.sX$counts)
   cdfY <- cumsum(df.sY$counts)
   ImodeX <- bX[df.sX$density==max(df.sX$density,na.rm=TRUE)]
   freqImodeX <- df.sX$counts[bX==ImodeX]
   ImodeY <- bY[df.sY$density==max(df.sY$density,na.rm=TRUE)]
   freqImodeY <- df.sY$counts[bY==ImodeY]
   mfun <- function(x,obj,g,cdf,Imode,freqImode,n=NULL) {
	if(x < Imode) return(1)
	if(is.null(n)) n <- length(obj$breaks)
	b <- obj$breaks
	id <- x > b[1:(n-1)] & x <= b[2:n]
	ind <- (1:(n-1))[id]
	return(1+g*round(cdf[ind]/freqImode, digits=0))
   } # end of internal 'mfun' function.
   N <- prod(xdim)
   mX <- mY <- numeric(N)+NA
   for(i in 1:N) {
	if(verbose & (i == 1 | i%%100==0)) cat(i, " ")
	if(tmpXc[i]>0) {
	   mX[i] <- mfun(x=tmpXc[i], obj=df.sX, g=gamma, cdf=cdfX, Imode=ImodeX, freqImode=freqImodeX,n=nX)
	   if(is.na(mX[i])) mX[i] <- 0
	}
	if(tmpYc[i]>0) {
	   mY[i] <- mfun(x=tmpYc[i], obj=df.sY, g=gamma, cdf=cdfY, Imode=ImodeY, freqImode=freqImodeY,n=nY)
	   if(is.na(mY[i])) mY[i] <- 0
	}
   } # end of for 'i' loop.
   if(verbose) cat("\n")
   mX <- mX[!is.na(mX)]
   mY <- mY[!is.na(mY)]
   tmp.sX <- cbind(rep(sX[,1], mX), rep(sX[,2], mX))
   sX <- tmp.sX
   tmp.sY <- cbind(rep(sY[,1], mY), rep(sY[,2], mY))
   sY <- tmp.sY
   # Ok, now try to fit the GMM using the EM algorithm.

   if(verbose) cat("\n", "Fitting GMM to verification/validation field.\n")
   fitX <- turboem(par=c(meanX,sdX,rep(1/K,K)), fixptfn = gmmEMstep, method="em", xy=sX, k=K)
   # fitX <- mvnormalmixEM(sX, lambda=rep(1/K,K), mu=meanX, sigma=sdX, k=K, verb=verbose)
   if(verbose) cat("\n", "Fitting GMM to forecast/model field.\n")
   # fitY <- mvnormalmixEM(sY, lambda=rep(1/K,K), mu=meanY, sigma=sdY, k=K, verb=verbose)
   fitY <- turboem(par=c(meanY,sdY,rep(1/K,K)), fixptfn = gmmEMstep, method="em", xy=sY, k=K)
   if(verbose) cat("\n", "Congratultions, the fits are made.\n")

    out <- list()
    attributes(out) <- atmp <- list(...)
    if(is.null(atmp$loc.byrow)) attr(out, "loc.byrow") <- FALSE

    out$fitX <- fitX
    out$fitY <- fitY
    out$initX <- c(meanX,sdX,rep(1/K,K))
    out$initY <- c(meanY,sdY,rep(1/K,K))
    out$sX <- sX
    out$sY <- sY
    out$k <- K
    out$Ax <- Ax
    out$Ay <- Ay
    out$xdim <- xdim

    class(out) <- "gmm2d"

    return(out)
} # end of 'gmm2d.default' function.

gmmEMstep <- function(p,xy,k) {
   np <- length(p)
   mu.x <- p[1:k]
   mu.y <- p[(k+1):(2*k)]
   sig <- matrix(p[(2*k+1):(6*k)], 4, k)
   lambdas <- p[(6*k+1):np]
   inverter <- function(x) {
     x <- matrix(x, 2, 2)
     x <- solve(x)
   return(c(x))
   } # end of internal-internal 'inverter' function.
   sig.inv <- apply(sig,2,inverter)
   sqrtdet <- function(x) {
     x <- matrix(x,2,2)
     x <- 1/(2*pi*sqrt(abs(det(x))))
     return(x)
   } # end of internal-internal 'sqrtdet' function.
   fac <- apply(sig,2,sqrtdet)
   nxy <- dim(xy)[1]
   f <- matrix(NA, nxy, k)
   effer <- function(z,a,S) a*exp(-((matrix(z,nrow=1) %*% S %*% matrix(z, ncol=1)))/2)
   for(i in 1:k) {
     z <- xy - cbind(rep(mu.x[i],nxy),rep(mu.y[i],nxy))
     SigInv <- matrix(sig.inv[,i], 2, 2)
     f[,i] <- apply(z, 1, effer, a=fac[i], S=SigInv)
   } # end of for 'i' loop.
   f <- f*matrix(lambdas,nxy,k,byrow=TRUE)
   Pxy <- rowSums(f, na.rm=TRUE)
   Pk  <- f/matrix(Pxy,nxy,k)
   PkSum <- colSums(Pk,na.rm=TRUE)
   # Now that we have the E step done, time to do the M step.
   mux.new <- colSums(Pk*matrix(xy[,1],nxy,k),na.rm=TRUE)/PkSum
   muy.new <- colSums(Pk*matrix(xy[,2],nxy,k),na.rm=TRUE)/PkSum
   sig11.new <- colSums(Pk*(matrix(xy[,1],nxy,k) - matrix(mu.x,nxy,k,byrow=TRUE))^2,na.rm=TRUE)/PkSum
   sig21.new <- colSums(Pk*(matrix(xy[,1],nxy,k) - matrix(mu.x,nxy,k,byrow=TRUE))*(matrix(xy[,2],nxy,k) - matrix(mu.y,nxy,k,byrow=TRUE)),na.rm=TRUE)/PkSum
   sig22.new <- colSums(Pk*(matrix(xy[,2],nxy,k) - matrix(mu.y,nxy,k,byrow=TRUE))^2,na.rm=TRUE)/PkSum
   lambdas.new <- colMeans(Pk,na.rm=TRUE)
   out <- c(mux.new, muy.new, c(rbind(sig11.new, sig21.new, sig21.new, sig22.new)), lambdas.new)
   return(out)
} # end of internal 'gmmEMstep' function.

gmmNegLogLik <- function(p,xy,k) {
   np <- length(p)
   mu.x <- p[1:k]
   mu.y <- p[(k+1):(2*k)]
   sig <- matrix(p[(2*k+1):(6*k)], 4, k)
   lambdas <- p[(6*k+1):np]
   inverter <- function(x) {
      x <- matrix(x, 2, 2)
      x <- solve(x)
      return(c(x))
   } # end of internal-internal 'inverter' function.
   sig.inv <- apply(sig,2,inverter)
   sqrtdet <- function(x) {
      x <- matrix(x,2,2)
      x <- 1/(2*pi*sqrt(abs(det(x))))
      return(x)
   } # end of internal-internal 'sqrtdet' function.
   fac <- apply(sig,2,sqrtdet)
   nxy <- dim(xy)[1]
   f <- matrix(NA, nxy, k)
   effer <- function(z,a,S) a*exp(-((matrix(z,nrow=1) %*% S %*% matrix(z, ncol=1)))/2)
   for(i in 1:k) {
      z <- xy - cbind(rep(mu.x[i],nxy),rep(mu.y[i],nxy))
      SigInv <- matrix(sig.inv[,i], 2, 2)
      f[,i] <- apply(z, 1, effer, a=fac[i], S=SigInv)
   } # end of for 'i' loop.
   f <- f*matrix(lambdas,nxy,k,byrow=TRUE)
   Pxy <- rowSums(f, na.rm=TRUE)
   return(-sum(log(Pxy)))
} # end of internal 'gmmNegLogLik' function.  

plot.gmm2d <- function(x, ..., col=c("gray", tim.colors(64)), zlim=c(0,1), horizontal=TRUE) {

    a <- attributes(x)
    loc.byrow <- a$loc.byrow

    if(is.null(a$projection)) proj <- FALSE
    else proj <- a$projection

    if(is.null(a$map)) domap <- FALSE
    else domap <- a$map

    xd <- x$xdim

    if(proj) {
	loc <- list(x=matrix(a$loc[,1], xd[1], xd[2], byrow=loc.byrow),
		    y=matrix(a$loc[,2], xd[1], xd[2], byrow=loc.byrow))
    }

    if(domap) {
	locr <- apply(a$loc, 2, range, finite=TRUE)
	ax <- list(x=pretty(round(a$loc[,1], digits=2)), y=pretty(round(a$loc[,2], digits=2)))
    }

    pX <- x$fitX$pars
    pY <- x$fitY$pars
    k  <- x$k

    if(!is.null(a$loc)) xloc <- a$loc
    else xloc <- cbind(rep(1:xd[1],xd[2]), rep(1:xd[2],each=xd[1]))

    tmp <- predict(object=x, x=xloc)
    predX <- matrix(x$Ax*tmp$predX, xd[1], xd[2])
    predY <- matrix(x$Ay*tmp$predY, xd[1], xd[2])
    # zl <- range(c(c(predX),c(predY)),finite=TRUE)

    if(is.null(a$msg)) par(mfrow=c(1,2))
    else par(mfrow=c(1,2), oma=c(0,0,2,0))

    if(is.null(a$data.name)) ptitle <- c("Verification", "Forecast")
    else {
	if(length(a$data.name) == 3) ptitle <- a$data.name[-1]
	else ptitle <- a$data.name
    }
     if(domap) {
 	if(proj) {
 	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
	    axis(2, at=ax$y, labels=ax$y)
 	    poly.image(loc$x, loc$y, predX, col=col, zlim=zlim, main=ptitle[1], add=TRUE)
 	    map(add=TRUE, lwd=1.5)
 	    map(add=TRUE, database="state")
 
 	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
	    axis(2, at=ax$y, labels=ax$y)
             poly.image(loc$x, loc$y, predY, col=col, zlim=zlim, main=ptitle[2], add=TRUE)
             map(add=TRUE, lwd=1.5)
             map(add=TRUE, database="state")
 	} else {
 	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
	    axis(2, at=ax$y, labels=ax$y)
 	    image(as.image(c(predX), x=a$loc, nrow=xd[1], ncol=xd[2], na.rm=TRUE), col=col, zlim=zlim, main=ptitle[1],
 		add=TRUE)
 	    map(add=TRUE, lwd=1.5)
             map(add=TRUE, database="state")
 
 	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
	    axis(2, at=ax$y, labels=ax$y)
 	    image(as.image(c(predY), x=a$loc, nrow=xd[1], ncol=xd[2], na.rm=TRUE), col=col, zlim=zlim, main=ptitle[2],
 		add=TRUE)
 	    map(add=TRUE, lwd=1.5)
             map(add=TRUE, database="state")
 	}
     } else {
	if(proj) {
    poly.image(loc$x, loc$y, predX, col=col, zlim=zlim, main=ptitle[1])
	    poly.image(loc$x, loc$y, predY, col=col, zlim=zlim, main=ptitle[2])
	} else {
           image(predX, col=col, zlim=zlim, main=ptitle[1])
           image(predY, col=col, zlim=zlim, main=ptitle[2])
	}
    }
    # image(as.image(c(predX), x=xloc, nx=x$xdim[1], ny=x$xdim[2], na.rm=TRUE), col=col, zlim=zlim, main=ptitle[1])
    # image(as.image(c(predY), x=xloc, nx=x$xdim[1], ny=x$xdim[2], na.rm=TRUE), col=col, zlim=zlim, main=ptitle[2])
    # image.plot(predX, legend.only=TRUE, horizontal=horizontal, col=col, zlim=zlim)

    if(!is.null(a$msg)) {
	title("")
	mtext(a$msg, line=0.05, outer=TRUE)
    }
    invisible()
} # end of 'plot.gmm2d' function.

predict.gmm2d <- function(object, ..., x) {

   pX <- object$fitX$pars
   pY <- object$fitY$pars

   np <- length(pX)
   k <- object$k

   # args <- list(...)

   # if(is.element("x",names(args))) x <- y <- args$x
   if(!missing(x)) y <- x
   else {

	x <- object$sX
	y <- object$sY

   }

   muXx <- pX[1:k]
   muXy <- pX[(k+1):(2*k)]
   muYx <- pY[1:k]
   muYy <- pY[(k+1):(2*k)]

   sigX <- matrix(pX[(2*k+1):(6*k)], 4, k)
   sigY <- matrix(pY[(2*k+1):(6*k)], 4, k)

   lambdasX <- pX[(6*k+1):np]
   lambdasY <- pY[(6*k+1):np]

   ifun <- function(x, mx, my, S, lam, k) {
	out <- numeric(k)+NA
	for(i in 1:k) {
	   z <- matrix(c(x[1] - mx[i], x[2] - my[i]), ncol=1)
	   Sig <- matrix(S[,i], 2, 2)
	   Sinv <- solve(Sig)
	   Sdet <- sqrt(abs(det(Sig)))
	   out[i] <- (lam[i]/(2*pi*Sdet))*exp(-(t(z) %*% Sinv %*% z)/2)
	} # end of for 'i' loop.
	return(sum(out,na.rm=TRUE))
   } # end of internal 'ifun' function.

   predX <- apply(x, 1, ifun, mx=muXx, my=muXy, S=sigX, lam=lambdasX, k=k)
   predY <- apply(y, 1, ifun, mx=muYx, my=muYy, S=sigY, lam=lambdasY, k=k)

   return(list(predX=predX, predY=predY))
} # end of 'predict.gmm2d' function.

summary.gmm2d <- function(object, ...) {
   out <- list()
   pX <- object$fitX$pars
   pY <- object$fitY$pars
   np <- length(pX)
   k <- object$k
   args <- list(...)
   if(is.element("silent",names(args))) silent <- args$silent
   else silent <- FALSE
   muXx <- pX[1:k]
   muXy <- pX[(k+1):(2*k)]
   muYx <- pY[1:k]
   muYy <- pY[(k+1):(2*k)]
   sigX <- matrix(pX[(2*k+1):(6*k)], 4, k)
   sigY <- matrix(pY[(2*k+1):(6*k)], 4, k)
   lambdasX <- pX[(6*k+1):np]
   lambdasY <- pY[(6*k+1):np]
   if(!silent) cat("\n", "Means for Verification field.\n")
   out$meanX <- tmp <- data.frame(mu.x=muXx, mu.y=muXy)
   if(!silent) print(tmp)
   if(!silent) cat("\n", "Means for Forecast field.\n")
   out$meanY <- tmp <- data.frame(mu.x=muYx, mu.y=muYy)
   if(!silent) print(tmp)
   if(!silent) cat("\n", "Covariances for Verification field.\n")
   out$covX <- list()
   out$covY <- list()
   for(i in 1:k) {
	out$covX[[i]] <- tmp <- matrix(sigX[,i],2,2)
	if(!silent) print(tmp)
   }
   if(!silent) cat("\n", "Covariances for Forecast field.\n")
   for(i in 1:k) {
	out$covY[[i]] <- tmp <- matrix(sigY[,i],2,2)
	if(!silent) print(tmp)
   }
   if(!silent) {
	cat("\n", "Mixture components (lambdas) for Verification field.\n")
        print(lambdasX)
        cat("\n", "Mixture components (lambdas) for Forecast field.\n")
        print(lambdasY)
   }
   out$lambdasX <- lambdasX
   out$lambdasY <- lambdasY
   # Now, actually make some calculations.
   efun <- function(x) return(eigen(x)$vectors[,1])
   eigX <- lapply(out$covX, efun)
   eigY <- lapply(out$covY, efun)
   e.tr <- matrix(NA, k, k)
   colnames(e.tr) <- as.character(1:k)
   rownames(e.tr) <- as.character(1:k)
   e.rot <- e.sc <- e.tr
   for(i in 1:k) for(j in 1:k) {
      e.tr[i,j] <- sqrt((muXx[i] - muYx[j])^2 + (muXy[i] - muYy[j])^2)
      e.rot[i,j] <- sum((180/pi)*acos(eigX[[i]]*eigY[[j]]))
      e.sc[i,j] <- object$Ay*lambdasY[j]/(object$Ax*lambdasX[i])
   } # end of double for 'i' and 'j' loop.
   if(!silent) {
	cat("\n", "Translation Errors between each component.\n")
	print(e.tr)
	cat("\n", "Rotation Errors between each component.\n")
        print(e.rot)
	cat("\n", "Scaling Errors between each component.\n")
        print(e.sc)
   }
   out$e.tr <- e.tr
   out$e.rot <- e.rot
   out$e.sc <- e.sc
   if(is.null(args$e1)) e1 <- 0.3
   else e1 <- args$e1
   if(is.null(args$e2)) e2 <- 100
   else e2 <- args$e2
   if(is.null(args$e3)) e3 <- 0.2
   else e3 <- args$e3
   if(is.null(args$e4)) e4 <- 90
   else e4 <- args$e4
   if(is.null(args$e5)) e5 <- 0.5
   else e4 <- args$e5
   e.overall <- e1*pmin(e.tr/e2,1) + e3*pmin(e.rot,180 - e.rot)/e4 + e5*(pmax(e.sc,1/e.sc)-1)
   out$e.overall <- e.overall
   if(!silent) {
	cat("\n", "Overall error between each component.\n")
	print(e.overall)
   }
   invisible(out)
} # end of 'summary.gmm2d' function.

print.gmm2d <- function(x, ...) {

    a <- attributes(x)

    cat("\n", "Gaussian Mixture Model Fits with ", x$k, " Models.\n")

    if(!is.null(a$data.name)) {

	print(a$data.name)
	dn <- a$data.name
	if(length(dn) == 3) dn <- dn[-1]
    } else dn <- c("verification", "forecast")

    if(!is.null(a$msg)) cat("\n", a$msg, "\n")

    y <- summary(x, silent=TRUE)

    cat(dn[1], " Means:\n") 
    print(y$meanX)

    cat("\n\n", dn[2], " Means:\n")
    print(y$meanY)

    invisible()
} # end of 'print.gmm2d' function.
