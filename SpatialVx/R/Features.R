summary.matched <- function(object, ...) {

    x <- object
    class( x ) <- "features"
    out <- summary( x )

    invisible( out )

} # end of 'summary.matched' function.

summary.features <- function(object,...) {

    args <- list(...)
    a <- attributes(object)

    if(is.null(args$silent)) silent <- FALSE
    else silent <- args$silent

    X <- object$X.feats
    Y <- object$Y.feats

    out <- list()
    b <- a
    b$names <- NULL
    attributes(out) <- b

    n <- length(X)
    m <- length(Y)

    if(!is.null(object$Xhat)) {
	Im1 <- object$X
	Im2 <- object$Xhat
	holdX <- matrix(NA, n, 7)
        holdY <- matrix(NA, m, 7)
        colnames(holdX) <- colnames(holdY) <- c("centroidX", "centroidY", "area", "OrientationAngle", "AspectRatio",
                                                "Intensity0.25", "Intensity0.9")
	wpr <- c("centroid", "area", "axis", "intensity")
	do.int <- TRUE
    } else {
	Im1 <- Im2 <- NULL
	holdX <- matrix(NA, n, 5)
        holdY <- matrix(NA, m, 5)
        colnames(holdX) <- colnames(holdY) <- c("centroidX", "centroidY", "area", "OrientationAngle", "AspectRatio")
	wpr <- c("centroid", "area", "axis")
	do.int <- FALSE
    }

    for(i in 1:n) {

      tmp <- FeatureProps(X[[i]], Im=Im1, which.props=wpr, loc = a$loc)
      holdX[i,1:2] <- c(tmp$centroid$x, tmp$centroid$y)
      holdX[i,3] <- tmp$area
      if(!is.null(c(tmp$axis$OrientationAngle$MajorAxis))) holdX[i,4] <- c(tmp$axis$OrientationAngle$MajorAxis)
      if(!is.null(tmp$axis$aspect.ratio)) holdX[i,5] <- tmp$axis$aspect.ratio
      if(do.int) holdX[i,6:7] <- c(tmp$intensity)

    }

    for(i in 1:m) {

      tmp <- FeatureProps(Y[[i]], Im=Im2, which.props=wpr, loc = a$loc)
      holdY[i,1:2] <- c(tmp$centroid$x, tmp$centroid$y)
      holdY[i,3] <- tmp$area
      if(!is.null(c(tmp$axis$OrientationAngle$MajorAxis))) holdY[i,4] <- c(tmp$axis$OrientationAngle$MajorAxis)
      if(!is.null(tmp$axis$aspect.ratio)) holdY[i,5] <- tmp$axis$aspect.ratio
      if(do.int) holdY[i,6:7] <- c(tmp$intensity)

    }

    if(!silent) {
	cat("\n", "Verification field (", object$data.name[1], ") feature properties:\n")
	print(holdX)
	cat("\n", "Forecast field (", object$data.name[2], ") feature properties:\n")
	print(holdY)
    }
    out$X <- holdX
    out$Y <- holdY
    class(out) <- "summary.features"

    invisible(out)

} # end of 'summary.features' function.

print.features <- function(x, ...) {

    a <- attributes(x)

    print( a$call )

    print(a$msg)
    print(a$data.name)

    if(length(a$data.name) == 3) {
	vxname <- a$data.name[2]
	fcstname <- a$data.name[3]
    } else {
	vxname <- a$data.name[1]
	fcstname <- a$data.name[2]
    }

    print(x$identifier.function)
    print(x$identifier.label)

    cat("\n", length(x$X.feats), " ", vxname, " features identified.\n")
    cat(length(x$Y.feats), " ", fcstname, " features identified.\n")

    if( !is.null( x$thresholds ) ) {

	cat("Thresholds used are:\n" )
	print( x$thresholds )

    }

    invisible()

} # end of 'print.features' function.

plot.features <- function(x, ..., type = c("both", "obs", "model")) {

    type <- tolower(type)
    type <- match.arg(type)

    if(type == "both") par(mfrow = c(1,2), oma = c(0,0,2,0))

    a <- attributes(x)
    loc.byrow <- a$loc.byrow

    if(length(a$data.name) == 3) {

        vxname <- a$data.name[2]
        fcstname <- a$data.name[3]

    } else {

        vxname <- a$data.name[1]
        fcstname <- a$data.name[2]

    }

    X <- x$X.labeled
    Y <- x$Y.labeled
    m <- max(c(c(X),c(Y)),na.rm=TRUE)

    icol <- c("white", rainbow(m))
    zl <- c(0,m)

    domap <- a$map
    proj <- a$projection
    xd <- a$xdim

    if(domap) {

	locr <- apply(a$loc, 2, range, finite=TRUE)
	ax <- list(x=pretty(round(a$loc[,1], digits=2)), y=pretty(round(a$loc[,2], digits=2)))

    }

    if(proj) loc <- list(x=matrix(a$loc[,1], xd[1], xd[2], byrow=loc.byrow),
		    y=matrix(a$loc[,2], xd[1], xd[2], byrow=loc.byrow))
    
    if(domap) {

	if(proj) {

	    if(is.element(type, c("both", "obs"))) {

	        map(xlim=locr[,1], ylim=locr[,2], type="n")
	        axis(1, at=ax$x, labels=ax$x)
	        axis(2, at=ax$y, labels=ax$y)

	        poly.image(loc$x, loc$y, X, col=icol, add=TRUE, xaxt="n", yaxt="n", zlim=zl)

	        map(add=TRUE, lwd=1.5)
	        map(add=TRUE, database="state")
	        title(vxname)

	    }

	    if(is.element(type, c("both", "model"))) {

	        map(xlim=locr[,1], ylim=locr[,2], type = "n")
	        axis(1, at=ax$x, labels=ax$x)
	        axis(2, at=ax$y, labels=ax$y)

                poly.image(loc$x, loc$y, Y, col = icol, add=TRUE, xaxt="n", yaxt="n", zlim=zl)

                map(add=TRUE, lwd=1.5)
                map(add=TRUE, database="state")
	        title(fcstname)

	    }

	} else {

	    if(is.element(type, c("both", "obs"))) {

	        map(xlim=locr[,1], ylim=locr[,2], type="n")
	        axis(1, at=ax$x, labels=ax$x)
	        axis(2, at=ax$y, labels=ax$y)

	        image(as.image(X, nx=xd[1], ny=xd[2], x=a$loc), col = icol, add=TRUE,
		    xaxt="n", yaxt="n", zlim=zl)

	        map(add=TRUE, lwd=1.5)
                map(add=TRUE, database="state")
	        title(vxname)

	    }

	    if(is.element(type, c("both", "model"))) {

	        map(xlim=locr[,1], ylim=locr[,2], type="n")
	        axis(1, at=ax$x, labels=ax$x)
	        axis(2, at=ax$y, labels=ax$y)

	        image(as.image(Y, nx=xd[1], ny=xd[2], x=a$loc), col = icol, add=TRUE,
		    xaxt="n", yaxt="n", zlim=zl)

	        map(add=TRUE, lwd=1.5)
                map(add=TRUE, database="state")

	        title(fcstname)

	    }

	}
    } else {

	if(proj) {

	    if(is.element(type, c("both", "obs"))) {

	        poly.image(loc$x, loc$y, X, col = icol, main=vxname, xaxt="n", yaxt="n", zlim=zl)

	    }

	    if(is.element(type, c("both", "model"))) {

	        poly.image(loc$x, loc$y, Y, col = icol, main=fcstname, xaxt="n", yaxt="n", zlim=zl)

	    }

	} else {

	    if(is.element(type, c("both", "obs"))) {

                image(X, col = icol, main=vxname, xaxt="n", yaxt="n", zlim=zl)

	    }

	    if(is.element(type, c("both", "model"))) {

                image(Y, col = icol, main=fcstname, xaxt="n", yaxt="n", zlim=zl)

	    }

	}
    }

    image.plot(Y, col = icol, zlim=zl, legend.only=TRUE, ...)

    if(type == "both" && !is.null(a$msg)) {

        title("")
        mtext(a$msg, line=0.05, outer=TRUE)

    }

    invisible()

} # end of 'plot.features' function.

plot.summary.features <- function(x, ...) {

   X <- x$X
   Y <- x$Y

   n <- dim(X)[1]
   m <- dim(Y)[1]
   k <- dim(X)[2]

   goodX <- sum(!is.na(X[,"OrientationAngle"]),na.rm=TRUE) > 2
   goodY <- sum(!is.na(Y[,"OrientationAngle"]),na.rm=TRUE) > 2

   if(k == 5) par(mfrow=c(2,2))
   else par(mfrow=c(4,2))
   
   hist(X[,"area"], col = "darkblue", breaks="FD", main="Verification \nFeature Area", xlab="area")
   hist(Y[,"area"], col = "darkorange", breaks="FD", main="Forecast \nFeature Area", xlab="area")

   if(goodX & goodY) {

	plot(X[,"OrientationAngle"], X[,"AspectRatio"], pch=19, col = "darkblue",
	    xlab="Major Axis Orientation Angle", ylab="Aspect Ratio")

   	points(Y[,"OrientationAngle"], Y[,"AspectRatio"], pch=19, col="darkorange")

   	legend("topright", legend=c("Verification", "Forecast"), pch=19, col=c("darkblue", "darkorange"), bty="n")
   }

   if(k==7) {
	hist(X[,"Intensity0.25"], col="darkblue", breaks="FD", main="Verification Feature \nIntensity (lower quartile)",
	    xlab="intensity")
	hist(X[,"Intensity0.9"], col="darkblue", breaks="FD", main="Verification Feature \nIntensity (0.9 quartile)",
	    xlab="intensity")
	hist(Y[,"Intensity0.25"], col="darkblue", breaks="FD", main="Forecast Feature \nIntensity (lower quartile)",
	    xlab="intensity")
        hist(Y[,"Intensity0.9"], col="darkblue", breaks="FD", main="Forecast Feature \nIntensity (0.9 quartile)",
	    xlab="intensity")
   }

    xl <- range(c(X[,"centroidX"],Y[,"centroidX"]), finite=TRUE)
    yl <- range(c(X[,"centroidY"],Y[,"centroidY"]), finite=TRUE)

    plot(X[,"centroidX"], X[,"centroidY"], pch=19, col="darkblue", xlim=xl, ylim=yl, xlab="", ylab="",
	main="Feature centroids")
    points(Y[,"centroidX"], Y[,"centroidY"], pch=19, col="darkorange")
    legend("topright", legend=c("Verification", "Forecast"), pch=19, col=c("darkblue", "darkorange"))

   invisible()
} # end of 'plot.summary.features' function.

disjointer <- function(x, method="C") {
   x[ x==0] <- NA
   if(any(!is.na(x))) {
	out <- as.im( x)
   	out <- connected(X=out, method=method)
   	out <- tiles(tess(image=out))
   } else out <- NULL
   return( out)
} # end of 'disjointer' function.

# threshfac <- function(object, fac=0.06666667, q=0.95, wash.out=NULL, thresh=NULL, idfun="disjointer",
#     time.point=1, model=1, ...) {
# 
#     theCall <- match.call()
# 
#     a <- attributes(object)
# 
#     ## Begin: Get the data sets
#     if(!missing(time.point) && !missing(model)) dat <- datagrabber(object, time.point=time.point, model=model)
#     else if(!missing(time.point)) dat <- datagrabber(object, time.point=time.point)
#     else if(!missing(model)) dat <- datagrabber(object, model=model)
#     else dat <- datagrabber(object)
# 
#     X <- dat$X
#     Y <- dat$Xhat
#     ## End: Get the data sets
# 
#     xdim <- a$xdim
#     Ix <- Iy <- matrix(0, xdim[1], xdim[2])
# 
#     if(is.null(thresh)) {
# 
# 	if(is.null(wash.out)){
# 
# 	   thresh <- quantile(c(X), probs=q)
# 	   thresh <- c(thresh, quantile(c(Y), probs=q))
# 
# 	} else {
# 
# 	   thresh <- quantile(c(X[X >= wash.out]), probs=q)
# 	   thresh <- c(thresh, quantile(c(Y[Y >= wash.out]), probs=q))
# 
# 	}
# 
# 	thresh <- thresh * fac
# 
#     } else if(length(thresh) == 1) thresh <- c(thresh, thresh)
# 
#     Ix[X >= thresh[1]] <- 1
#     Iy[Y >= thresh[2]] <- 1
# 
#     X.feats <- do.call(idfun, c(list(x=Ix), list(...)))
#     Y.feats <- do.call(idfun, c(list(x=Iy), list(...)))
# 
#     Xlab <- Ylab <- matrix(0, xdim[1], xdim[2])
# 
#     if(!is.null(X.feats)) for( i in 1:length( X.feats)) Xlab[X.feats[[i]][["m"]]] <- i
#     else X.feats <- NULL
# 
#     if(!is.null(Y.feats)) for( j in 1:length( Y.feats)) Ylab[Y.feats[[j]][["m"]]] <- j
#     else Y.feats <- NULL
# 
#     # if(is.null(X.feats)) warning("threshfac: No values above threshold in verification field.")
#     # if(is.null(Y.feats)) warning("threshfac: No values above threshold in forecast field.")
# 
#     out <- list()
#     attributes(out) <- a
# 
#     out$X <- X
#     out$Xhat <- Y
#     out$X.feats <- X.feats
#     out$Y.feats <- Y.feats
#     out$X.labeled <- Xlab
#     out$Y.labeled <- Ylab
#     out$identifier.function <- "threshfac"
#     out$identifier.label <- "Threshold"
#     names( thresh ) <- c( "Observed", "Forecast" )
#     out$threshold <- thresh
# 
#     attr(out, "time.point") <- time.point
#     attr(out, "model") <- model
# 
#     if(length(a$data.name) == a$nforecast + 2) {
# 
#         dn <- a$data.name[-(1:2)]
#         vxname <- a$data.name[1:2]
# 
#     } else {
# 
#         dn <- a$data.name[-1]
#         vxname <- a$data.name[1]
# 
#     }
# 
#     if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
#     else model.num <- model
# 
#     attr(out, "data.name") <- c(vxname, dn[model.num])
# 
#     attr( out, "call") <- theCall
# 
#     class(out) <- "features"
#     return(out)
# 
# } # end of 'threshfac' function.

saller <- function(x, d=NULL, distfun = "rdist", ...) {

    out <- list()
    a <- attributes(x)
    if(!is.null(a$names)) a$names <- NULL
    attributes(out) <- a

    tmp <- x

    y <- tmp$Y.feats
    x <- tmp$X.feats

    binX <- im(tmp$X.labeled)
    binX <- solutionset(binX > 0)
    binY <- im(tmp$Y.labeled)
    binY <- solutionset(binY > 0)

    # The verification set (images).
    X <- tmp$X
    Y <- tmp$Xhat

    xdim <- dim(X)

    # Amplitude
    DomRmod <- mean(Y,na.rm=TRUE)
    DomRobs <- mean(X,na.rm=TRUE)
    A <- 2*(DomRmod - DomRobs)/(DomRmod + DomRobs)
    out$A <- A

    # Location
    if(is.null(d)) d <- max(xdim, na.rm=TRUE)
    num <- centdist(binY,binX, distfun = distfun, loc = a$loc, ...)
    L1 <- num/d
    intRamt <- function(id,x) return(sum(x[id$m],na.rm=TRUE))
    RnMod <- as.numeric(unlist(lapply(y, intRamt, x=Y)))
    RnObs <- as.numeric(unlist(lapply(x, intRamt, x=X)))
    xRmodN <- as.numeric(unlist(lapply(y, centdist, y=binY)))
    xRobsN <- as.numeric(unlist(lapply(x, centdist, y=binX)))
    RmodSum <- sum( RnMod, na.rm=TRUE)
    RobsSum <- sum( RnObs, na.rm=TRUE)
    rmod <- sum( RnMod*xRmodN, na.rm=TRUE)/RmodSum
    robs <- sum( RnObs*xRobsN, na.rm=TRUE)/RobsSum
    L2 <- 2*abs(rmod - robs)/d
    out$L1 <- L1
    out$L2 <- L2
    out$L <- L1 + L2

    # Structure
    Rmaxer <- function(id, x) return(max(x[id$m], na.rm=TRUE))
    RnMaxMod <- as.numeric(unlist(lapply(y, Rmaxer, x=Y)))
    RnMaxObs <- as.numeric(unlist(lapply(x, Rmaxer, x=X)))
    VmodN <- RnMod/RnMaxMod
    VobsN <- RnObs/RnMaxObs
    Vmod <- sum(RnMod*VmodN, na.rm=TRUE)/RmodSum
    Vobs <- sum(RnObs*VobsN, na.rm=TRUE)/RobsSum
    out$S <- 2*(Vmod - Vobs)/(Vmod + Vobs)

    class(out) <- "saller"
    return(out)

} # end of 'saller' function.

print.saller <- function(x, ...) {

    a <- attributes(x)
    cat(a$msg, "\n")
    print(a$data.name)
    cat("\n\n")
    b <- c(x$S, x$A, x$L)
    names(b) <- c("S", "A", "L")
    print(b)

    invisible(b)

} # end of 'print.saller' function.

summary.saller <- function(object,...) {

   # args <- list(...)
    a <- attributes(object)
    print(a$msg)
    print(a$data.name)

    cat("\n", "Structure Component (S): ", object$S, "\n")
    cat("\n", "Amplitude Component (A): ", object$A, "\n")
    cat("\n", "Location Component (L): ", object$L, "\n")

    invisible()

} # end of 'summary.saller' function.

plot.saller <- function(x, ...) {

   invisible() 

} # end of 'plot.saller' function.

centdist <- function(x,y, distfun = "rdist", loc = NULL, ...) {

    # xcen <- centroid.owin(x)
    # ycen <- centroid.owin(y)

    xcen <- FeatureProps(x = x, which.props = "centroid", loc = loc)
    ycen <- FeatureProps(x = y, which.props = "centroid", loc = loc)

    x1 <- matrix(c(xcen$centroid$x, xcen$centroid$y), 1, 2)
    x2 <- matrix(c(ycen$centroid$x, ycen$centroid$y), 1, 2)

    out <- c(do.call(distfun, c(list(x1 = x1, x2 = x2), list(...))))

    # return(sqrt((xcen$x - ycen$x)^2 + (xcen$y - ycen$y)^2))
    return(out)

} # end of 'centdist' function.

plot.matched <- function(x, ..., type = c("both", "obs", "model")) {

    a <- attributes(x)
    loc.byrow <- a$loc.byrow

    args <- list(...)

    type <- tolower(type)
    type <- match.arg(type)

    matches <- x$matches
    if(!is.null(x$implicit.merges)) mer <- x$implicit.merges
    else mer <- x$merges

    xdim <- dim(x$X.labeled)

    # Need to fill in values for X and Xhat
    # where integers correspond to the first n
    # matched objects and n + 1 for all unmatched
    # objects.

    if(any(dim(matches) == 0)) {

	n <- 0

	X <- x$X.labeled
	X[X > 0] <- 1

	Xhat <- x$Y.labeled
	Xhat[Xhat > 0] <- 1

	icol <- c("white", "gray")

    } else {

	X <- Xhat <- matrix(0, xdim[1], xdim[2])

        if(is.null(mer)) { 

	    n <- dim(matches)[1]

	    for(i in 1:n) {

		k <- matches[i,"Observed"]
		j <- matches[i, "Forecast"]

		look <- x$X.feats[[ k ]]
		look <- look$m

		X[look] <- i

		look <- x$Y.feats[[ j ]]
		look <- look$m

		Xhat[look] <- i

	    } # end of for 'i' loop.

        } else {

	    n <- length(mer)

	    for(i in 1:n) {

		mi <- mer[[i]] # mi is like matches.

		for(jj in 1:dim(mi)[1]) {

		    k <- mi[jj, "Observed"]
		    j <- mi[jj, "Forecast"]

		    look <- x$X.feats[[ k ]]
		    look <- look$m

		    X[look] <- i

		    look <- x$Y.feats[[ j ]]
		    look <- look$m

		    Xhat[look] <- i

		}

	    } # end of for 'i' loop.

        } # end of if else 'implicit.merges/merges' stmts.

	oun <- x$unmatched$X
	unolen <- length(oun)
	fcun <- x$unmatched$Xhat
	unflen <- length(fcun)

	if(unolen > 0) {

	    for(i in 1:unolen) {

		look <- x$X.feats[[ oun[i] ]]
		look <- look$m

		X[look] <- n + 1

	    } # end of for 'i' loop.

	} # end of if any unmatched observed features stmt.

	if(unflen > 0) {

	    for(i in 1:unflen) {

		look <- x$Y.feats[[ fcun[i] ]]
		look <- look$m

		Xhat[look] <- n + 1

	    } # end of for 'i' loop.

	} # end of if any unmatched forecast features stmt.

	icol <- c("white", rainbow(n), "gray")

    } # end of if no matches stmt.

    if(!is.null(a$data.name)) {

        dn <- a$data.name

        if(length(dn) == 3) {

            vxname <- dn[2]
            fcstname <- dn[3]

        } else {

            vxname <- dn[1]
            fcstname <- dn[2]

        }

        X.name <- paste(vxname, "\nFeature Field", sep="")
        Xhat.name <- paste(fcstname, "\nFeature Field", sep="")

    } else {

            X.name <- "Verification\nFeature Field"
            Xhat.name <- "Forecast\nFeature Field"

    } # end of if '!is.null(a$data.name)' stmts.

    if(type == "both") par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))

    if (is.null(a$projection)) proj <- FALSE
    else proj <- a$projection

    if (is.null(a$map)) domap <- FALSE
    else domap <- a$map

    if(proj) loc <- list(x = matrix(a$loc[, 1], xdim[1], xdim[2], byrow = loc.byrow), 
            		y = matrix(a$loc[, 2], xdim[1], xdim[2], byrow = loc.byrow))

    zl <- c(0, n + 1)

    if(domap) {

        locr <- apply(a$loc, 2, range, finite = TRUE)

        ax <- list(x = pretty(round(a$loc[, 1], digits = 2)), 
            	   y = pretty(round(a$loc[, 2], digits = 2)))

	if(proj) {

	   if(is.element(type, c("both", "obs"))) {

	        map(xlim = locr[, 1], ylim = locr[, 2], type = "n")

                axis(1, at = ax$x, labels = ax$x)
                axis(2, at = ax$y, labels = ax$y)

                poly.image(loc$x, loc$y, X, add = TRUE, col = icol, zlim = zl)
                map(add = TRUE, lwd = 1.5)
                map(add = TRUE, database = "state")
                title(X.name)

	    }

	    if(is.element(type, c("both", "model"))) {

                map(xlim = locr[, 1], ylim = locr[, 2], type = "n")
                axis(1, at = ax$x, labels = ax$x)
                axis(2, at = ax$y, labels = ax$y)

                poly.image(loc$x, loc$y, Xhat, add = TRUE, col = icol, zlim = zl)
                map(add = TRUE, lwd = 1.5)
                map(add = TRUE, database = "state")
                title(Xhat.name) 

	    }

	} else {

	    if(is.element(type, c("both", "obs"))) {

	        map(xlim = locr[, 1], ylim = locr[, 2], type = "n")
                axis(1, at = ax$x, labels = ax$x)
                axis(2, at = ax$y, labels = ax$y)

                image(as.image(X, nx = xdim[1], ny = xdim[2], x = a$loc, na.rm = TRUE),
		    col = icol, zlim = zl, add = TRUE)
                map(add = TRUE, lwd = 1.5)
                map(add = TRUE, database = "state")
                title(X.name)

	    }

	    if(is.element(type, c("both", "model"))) {

                map(xlim = locr[, 1], ylim = locr[, 2], type = "n")
                axis(1, at = ax$x, labels = ax$x)
                axis(2, at = ax$y, labels = ax$y)

                image(as.image(Xhat, nx = xdim[1], ny = xdim[2], x = a$loc, na.rm = TRUE),
		    col = icol, zlim = zl, add = TRUE)
                map(add = TRUE, lwd = 1.5)
                map(add = TRUE, database = "state")
                title(Xhat.name)

	    }

	} # end of if else 'proj' stmt.

    } else {

	if (proj) {

	    if(is.element(type, c("both", "obs"))) {

                poly.image(loc$x, loc$y, X, add = TRUE, col = icol, zlim = zl)
                title(X.name)

	    }

	    if(is.element(type, c("both", "model"))) {
 
                poly.image(loc$x, loc$y, Xhat, add = TRUE, col = icol, zlim = zl)
                title(Xhat.name)

	    }

        } else {

	    if(is.element(type, c("both", "obs"))) {

                image(X, col = icol, zlim = zl, main = X.name)

	    }

	    if(is.element(type, c("both", "model"))) {

                image(Xhat, col = icol, zlim = zl, main = Xhat.name)

	    }

        } # end of if else 'proj' stmt.


    } # end of if else 'domap' stmts.

    image.plot(X, col = icol, zlim = zl, legend.only = TRUE, ...)

    if(!is.null(a$msg)) {

	title("")
	mtext(a$msg, line = 0.05, outer = TRUE)

    }

    invisible()

} # end of 'plot.matched' function.

print.matched <- function(x, ...) {

    a <- attributes(x)
    print(x$match.message)

    print(x$match.type)

    if(!is.null(x$criteria)) {
	
	if(x$criteria == 1) cat("Distance criteria based on sum of the object sizes.\n")
	else if(x$criteria == 2) cat("Distance criteria based on the average of the object sizes.\n")
	else if(x$criteria == 3) cat("Distance criteria based on a constant threshold given by ", x$const, "\n")

    }

    cat("Matched objects.\n")
    print(x$matches)

    cat("Unmatched Objects.\n")
    print(x$unmatched)

    if(!is.null(x$implicit.merges)) {

	cat("Objects are not merged, but implicitly defined merges would be:\n")

	for(i in 1:length(x$implicit.merges)) {

	    cat("New object ", i, ":\n")
	    print(x$implicit.merges[[i]])

	}

    }

    invisible()

} # end of 'print.matched' function.

centmatch <- function(x, criteria = 1, const = 14, distfun = "rdist", areafac = 1,
    verbose = FALSE, ...) {

    if(class(x) != "features") stop("centmatch: invalid object, x or y type.")

    out <- x

    out$match.message <- "Matching based on centroid distances using centmatch function."
    out$match.type <- "centmatch"

    out$criteria <- criteria
    if(criteria == 3) out$const <- const
    else out$const <- NULL

    a <- attributes(x)

    if(distfun == "rdist.earth") {

	loc <- a$loc
	if(is.null(loc)) warning("Using rdist.earth, but lon/lat coords are not available. Can pass them as an attribute to x called loc.")

    } else loc <- NULL

    xdim <- dim(x$X.labeled)

    # Get the list of "owin" class objects defining individual features
    # for each field.
    Y <- x$Y.feats
    X <- x$X.feats


    m <- length(Y)
    n <- length(X)

    if(criteria != 3) {

      Ax <- numeric(n)
      Ay <- numeric(m)

    }

    # matrix to hold the centroid distances between each pair of m model and n observed objects.
    Dcent  <- matrix(NA,m,n)

    if(verbose) {

	if(criteria != 3) cat("\n", "Looping through each feature in each field to find the centroid differences.\n")
	else cat("\n", "Looping through each feature in each field to find the areas and centroid differences.\n")

    }

    for(i in 1:m) {

	if(verbose) cat(i, "\n")

	if(criteria != 3) {

	   tmpy <- FeatureProps(x=Y[[i]], which.props=c("centroid", "area"), areafac = areafac, loc = loc)
	   Ay[i] <- sqrt(tmpy$area)

	} else tmpy <- FeatureProps(x=Y[[i]], which.props="centroid", areafac = areafac, loc = loc)

	ycen <- tmpy$centroid

	for(j in 1:n) {

	   if(verbose) cat(j, " ")
	   if(criteria != 3) {

		tmpx <- FeatureProps(x=X[[j]], which.props=c("centroid", "area"), areafac = areafac, loc = loc)
	   	Ax[j] <- sqrt(tmpx$area)

	   } else tmpx <- FeatureProps(x=X[[j]], which.props="centroid", areafac = areafac, loc = loc)

	   xcen <- tmpx$centroid
	   # Dcent[i,j] <- sqrt((xcen$x - ycen$x)^2 + (xcen$y - ycen$y)^2)
	   Dcent[i, j] <- do.call(distfun, c(list(x1 = matrix(c(xcen$x, xcen$y), 1, 2),
					    x2 = matrix(c(ycen$x, ycen$y), 1, 2), ...)))

	} # end of for 'j' loop.

	if(verbose) cat("\n")

    } # end of for 'i' loop.

    if(criteria != 3) {

	Ay <- matrix( rep(Ay,n), m, n)
   	Ax <- matrix( rep(Ax,m), m, n, byrow=TRUE)

    }

    if(criteria == 1) Dcomp <- Ay + Ax
    else if(criteria == 2) Dcomp <- (Ax + Ay)/2
    else if(criteria == 3) Dcomp <- matrix(const,m,n)
    else stop("centmatch: criteria must be 1, 2 or 3.")

    DcompID <- Dcent < Dcomp

    any.matched <- any(DcompID)

    FobjID <- matrix( rep(1:m, n), m, n)
    OobjID <- t(matrix( rep(1:n, m), n, m))

    fmatches <- cbind(c(FobjID)[DcompID], c(OobjID)[DcompID])
    colnames(fmatches) <- c("Forecast", "Observed")

    if(dim(fmatches)[ 1 ] > 1) {

        # Check for multiple object pairs.
        pcheck <- paste(fmatches[,1], fmatches[,2], sep="-")
        dupID <- duplicated(pcheck)
        if(any(dupID)) fmatches <- fmatches[!dupID, , drop = FALSE]

        # Now, put the objects in order according to the forecast objects.
        oID <- order(fmatches[,1])
        fmatches <- fmatches[oID, , drop = FALSE]

    }

    if(is.null(dim(fmatches)) && length(fmatches) == 2) fmatches <- matrix(fmatches, ncol = 2)
    out$matches <- fmatches

    if(any.matched) {

	funmatched <- (1:m)[!is.element(1:m, fmatches[,1])]
	vxunmatched <- (1:n)[!is.element(1:n, fmatches[,2])]

	# Find implicit merges.  That is, objects from one field
	# can be matched multiple times to objects in the other field.
	# While these are not to be considered merged, it might make
	# sense to consider them merged.  Defining those merges here
	# should make things easier later (e.g., when plotting) uses.

	matchlen <- dim(fmatches)[1]
	fuq <- unique(fmatches[,1])
	flen <- length(fuq)
	ouq <- unique(fmatches[,2])
	olen <- length(ouq)

	if(matchlen == flen && matchlen > olen) {

	    if(verbose) cat("Multiple observed features are matched to one or more forecast feature(s).  Determining implicit merges.\n")

	} else if(matchlen > flen && matchlen == olen) {

	    if(verbose) cat("Multiple forecast features are matched to one or more observed feature(s).  Determining implicit merges.\n")

	} else if(matchlen > flen && matchlen > olen) {

	    if(verbose) cat("Multiple matches have been found between features in each field.  Determining implicit merges.\n")

	} else if(matchlen == flen && matchlen == olen) {

	    if(verbose) cat("No multiple matches were found.  Thus, no implicit merges need be considered.\n")

	} # end of if else which fields have multiple matches stmts.

	implicit.merges <- MergeIdentifier(fmatches)

    } else {

	if(verbose) cat("No objects matched.\n")

	implicit.merges <- NULL

	funmatched <- 1:m
	vxunmatched <- 1:n

    }

    out$unmatched <- list(X = vxunmatched, Xhat = funmatched)

    out$implicit.merges <- implicit.merges

    out$criteria.values <- Dcomp
    out$centroid.distances <- Dcent

    class(out) <- "matched"
    return(out)

} # end of 'centmatch' function.

MergeIdentifier <- function(x) {
    
    if(any(dim(x) == 0)) return(NULL)

    matchlen <- dim(x)[1]
    
    fuq <- unique(x[,1])
    ouq <- unique(x[,2])

    flen <- length(fuq)
    olen <- length(ouq)

    out <- list()

    if(matchlen == flen && matchlen > olen) for(i in 1:olen) out[[ i ]] <- x[ x[,2] == ouq[ i ], , drop = FALSE]
    else if(matchlen > flen && matchlen == olen) for(i in 1:flen) out[[ i ]] <- x[ x[,1] == fuq[ i ], , drop = FALSE]
    else if(matchlen > flen && matchlen > olen) {

        if(matchlen > 1) {

            idx <- integer(matchlen) + NA
            idx[1] <- 1

            for(i in 2:matchlen) {

                idF <- x[1:(i - 1), 1] == x[i, 1]
                idO <- x[1:(i - 1), 2] == x[i, 2]

                if(any(idF)) idx[i] <- idx[ (1:(i - 1))[idF] ][1]
                else if(any(idO)) idx[i] <- idx[ (1:(i - 1))[idO] ][1]
                else idx[i] <- max(idx[1:(i - 1)], na.rm = TRUE) + 1

            } # end of for 'i' loop.

        } else return(NULL) # just in case...

        for(j in 1:max(idx, na.rm = TRUE)) out[[ j ]] <- x[ idx == j, , drop = FALSE]
    
    } else if(matchlen == flen && matchlen == olen) return(NULL)

    return(out)

} # end of 'MergeIdentifier' function.

FeatureProps <- function(x, Im = NULL, which.props = c("centroid", "area", "axis", "intensity"),
    areafac = 1, q = c(0.25, 0.90), loc = NULL, ...) {

    out <- list()

    if(is.element("centroid", which.props)) {

        if(is.null(loc)) {

            xd <- dim(x$m)
            loc <- cbind(rep(1:xd[1], xd[2]), rep(1:xd[2], each = xd[1]))

        }

        xcen <- mean(loc[c(x$m), 1], na.rm = TRUE)
        ycen <- mean(loc[c(x$m), 2], na.rm = TRUE)
        centroid <- list(x = xcen, y = ycen)
        # out$centroid <- centroid.owin(x)

        out$centroid <- centroid

    } # end of if find centroid stmt.

    if(is.element("area", which.props)) out$area <- sum(colSums(x$m, na.rm=TRUE), na.rm=TRUE)*areafac

    if(is.element("axis", which.props)) out$axis <- FeatureAxis(x=x,fac=areafac, ...)

    if(is.element("intensity", which.props)) {

	ivec <- matrix(NA, ncol=length(q), nrow=1)
	colnames(ivec) <- as.character(q)
	ivec[1,] <- quantile(c(Im[x$m]), probs=q)
	out$intensity <- ivec

    }

    return(out)

} # end of 'FeatureProps' function.

FeatureComps <- function(Y, X, which.comps=c("cent.dist", "angle.diff", "area.ratio",
					     "int.area", "bdelta", "haus", "ph", 
					     "med", "msd", "fom", "minsep", "bearing"),
		sizefac=1, alpha=0.1, k=4, p=2, c=Inf, distfun="distmapfun", deg=TRUE, aty="compass", loc = NULL, ...) {

   id1 <- is.element(c("cent.dist", "angle.diff", "area.ratio", "int.area", "bearing"), which.comps)
   if(any(id1)) {
	list1 <- character(0)
	if( any(is.element(c("cent.dist","bearing"), which.comps))) list1 <- c(list1, "centroid")
	if( any(is.element(c("area.ratio","int.area"), which.comps))) list1 <- c(list1, "area")
	if( is.element("angle.diff", which.comps)) list1 <- c(list1, "axis")
   }
   id2 <- is.element(c("ph", "med", "msd", "fom", "minsep"), which.comps)
   if(any(id2)) list2 <- c("ph", "med", "msd", "fom", "minsep")[id2]

   if(any(id1)) {

	Xsingle.props <- FeatureProps(x=X, which.props=list1, areafac=sizefac^2, loc = loc)
	Ysingle.props <- FeatureProps(x=Y, which.props=list1, areafac=sizefac^2, loc = loc)

   }

   if(any(id2)) {

	# Xim <- as.im(X)
        # Xim <- solutionset(X == 1)

	# Yim <- as.im(Y)
        # Yim <- solutionset(Y == 1)

	out <- locperf(X=X, Y=Y, which.stats=list2, alpha=alpha, k=k, distfun=distfun, ...)

   } else out <- list()

   if(is.element("cent.dist", which.comps)) {

	Xcent.x <- Xsingle.props$centroid$x
	Xcent.y <- Xsingle.props$centroid$y
	Ycent.x <- Ysingle.props$centroid$x
	Ycent.y <- Ysingle.props$centroid$y
	out$cent.dist <- sqrt( (Ycent.x - Xcent.x)^2 + (Ycent.y - Xcent.y)^2)*sizefac

   } # end of if do centroid distance stmt.

   if(is.element("angle.diff", which.comps)) {

	phiX <- Xsingle.props$axis$OrientationAngle$MajorAxis*pi/180
	phiY <- Ysingle.props$axis$OrientationAngle$MajorAxis*pi/180
	out$angle.diff <- abs(atan2(sin(phiX-phiY),cos(phiX-phiY))*180/pi)

   } # end of if do angle difference stmts.

   if(any(is.element(c("area.ratio","int.area"), which.comps))) {

	Xa <- Xsingle.props$area
        Ya <- Ysingle.props$area

   } # end of if do any area calculations stmt.

   if(is.element("area.ratio", which.comps)) out$area.ratio <- min(Xa,Ya)/max(Xa,Ya)

   if(is.element("int.area", which.comps)) {

	denom <- (Xa + Ya)/2
	XY <- intersect.owin(X,Y)
	XYa <- FeatureProps(XY, which.props="area", areafac=sizefac^2, loc = loc)$area
	out$int.area <- XYa/denom

   } # end of if do area ratio stmt.

   if(is.element("bearing", which.comps)) out$bearing <- bearing(cbind(Ysingle.props$centroid$x,Ysingle.props$centroid$y),
									cbind(Xsingle.props$centroid$x,Xsingle.props$centroid$y),
								deg=deg, aty=aty)

   if(is.element("bdelta", which.comps)) out$bdelta <- deltametric(X,Y, p=p, c=c)
   if(is.element("haus", which.comps)) out$haus <- deltametric(X,Y,p=Inf,c=Inf)

   class(out) <- "FeatureComps"
   return(out)

} # end of 'FeatureComps' function.

FeatureAxis <- function(x, fac=1, flipit=FALSE, twixt=FALSE) {

   out <- list()
   # out$feature.name <- as.character(substitute(x))
   if( flipit) x <- flipxy(x)
   out$z <- x
   ch <- convexhull(x)
   out$chull <- ch
   pts <- unname(cbind(ch$bdry[[1]][["x"]], ch$bdry[[1]][["y"]]))
   out$pts <- pts
   axfit <- sma(y~x, data.frame(x = pts[,1], y = pts[,2]))
   axis.x <- c(axfit$from[[1]], axfit$to[[1]])
   a <- axfit$coef[[1]][1,1]
   b <- axfit$coef[[1]][2,1]
   axis.y <- a + b * axis.x

   if(any(c(is.na(axis.x),is.na(axis.y)))) return(NULL)
   axwin <- owin(xrange=range(axis.x), yrange=range(axis.y))
   MajorAxis <- as.psp(data.frame(x0 = axis.x[1], y0 = axis.y[1], x1 = axis.x[2], y1 = axis.y[2]), window = axwin)

   theta <- angles.psp(MajorAxis)
   if((0 <= theta) & (theta <= pi/2)) theta2 <- pi/2 - theta
   else theta2 <- 3 * pi/2 - theta
   tmp <- rotate(ch,theta2)
   tmp <- boundingbox(tmp)
   l <- tmp$xrange[2] - tmp$xrange[1]
   theta <- theta * 180/pi

   if(twixt) {

	if((theta > 90) & (theta <= 270)) theta <- theta - 180
   	else if((theta > 270) & (theta <= 360)) theta <- theta - 360

   } # end of if force theta to be between +/- 90 degrees.

   MidPoint <- midpoints.psp(MajorAxis)

   r <- lengths.psp(MajorAxis) * fac

   phi <- angles.psp(rotate(MajorAxis, pi/2))

   MinorAxis <- as.psp(data.frame(xmid = MidPoint$x, ymid = MidPoint$y, length = l / fac, angle = phi), window = axwin)
   phi <- phi * 180/pi

   out$MajorAxis <- MajorAxis
   out$MinorAxis <- MinorAxis
   out$OrientationAngle <- list(MajorAxis=theta, MinorAxis=phi)
   out$aspect.ratio <- l/r
   out$MidPoint <- MidPoint
   out$lengths <- list(MajorAxis=r, MinorAxis=l)
   out$sma.fit <- axfit
   class(out) <- "FeatureAxis"
   return(out)

} # end of 'FeatureAxis' function.

plot.FeatureAxis <- function(x, ..., zoom = FALSE) {
   args <- list(...)
   if(!zoom) plot( x$z, col="darkblue", main="", ...)
   else {

	z <- x$z
	bb <- boundingbox(z)
	z <- rebound(z, rect = bb)
	plot(z, col="darkblue", main="", ...)

   }

   plot( x$chull, add=TRUE)
   plot( x$MajorAxis, add=TRUE, col="yellow", lwd=1.5)
   plot( x$MajorAxis, add=TRUE, lty=2, lwd=1.5)
   plot( x$MinorAxis, add=TRUE, col="yellow", lwd=1.5)
   plot( x$MinorAxis, add=TRUE, lty=2, lwd=1.5)
   plot( x$MidPoint, add=TRUE, col="darkorange")

   invisible()

} # end of 'plot.FeatureAxis' function.

summary.FeatureAxis <- function(object, ...) {
   cat("\n", "Mid-point of Axis is at:\n")
   print(paste("(", object$MidPoint$x, ", ", object$MidPoint$y, ")", sep=""))
   cat("\n", "Major Axis length = ", object$lengths$MajorAxis, "\n")
   cat("\n", "Major Axis Angle = ", object$OrientationAngle$MajorAxis, " degrees\n")
   cat("\n", "Minor Axis length = ", object$lengths$MinorAxis, "\n")
   cat("\n", "Minor Axis Angle = ", object$OrientationAngle$MinorAxis, " degrees\n")
   cat("\n", "Aspect ratio = ", object$aspect.ratio, "\n")
   cat("\n\n", "sma fit summary (see help file for function sma from package smatr)\n\n")
   print(summary(object$sma.fit))
   invisible()
} # end of 'summary.FeatureAxis' function.

FeatureMatchAnalyzer <- function(x, which.comps=c("cent.dist", "angle.diff", "area.ratio", "int.area",
                    "bdelta", "haus", "ph", "med", "msd", "fom", "minsep", "bearing"), 
		    sizefac=1, alpha=0.1, k=4, p=2, c=Inf, distfun="distmapfun", ...) {

    UseMethod("FeatureMatchAnalyzer", x)

} # end of 'FeatureMatchAnalyzer' function.

FeatureMatchAnalyzer.matched <- function(x, which.comps=c("cent.dist", "angle.diff", "area.ratio", "int.area",
                    "bdelta", "haus", "ph", "med", "msd", "fom", "minsep", "bearing"),
                    sizefac=1, alpha=0.1, k=4, p=2, c=Inf, distfun="distmapfun", ...) {

    class(x) <- paste(class(x), x$match.type, sep = ".")
    UseMethod("FeatureMatchAnalyzer", x)

} # end of 'FeatureMatchAnalyzer.matched' function.


FeatureMatchAnalyzer.matched.centmatch <- function(x, which.comps=c("cent.dist", "angle.diff", "area.ratio", "int.area",
                    "bdelta", "haus", "ph", "med", "msd", "fom", "minsep", "bearing"), 
                    sizefac=1, alpha=0.1, k=4, p=2, c=Inf, distfun="distmapfun", ...) {

    a <- attributes(x)
    a$names <- NULL

    n <- dim(x$matches)[1]

    if(n > 0) {
	
        out <- list()
	attributes(out) <- a

	loc <- a$loc

	Xfeats <- x$X.feats
	Yfeats <- x$Y.feats

	for(i in 1:n) {

	    j <- x$matches[i, "Forecast"]
	    k <- x$matches[i, "Observed"]
	    out[[i]] <- FeatureComps(Y=Yfeats[[ j ]], X=Xfeats[[ k ]], which.comps=which.comps,
                                        sizefac=sizefac, alpha=alpha, k=k, p=p, c=c, distfun=distfun, loc = loc, ...)

	} # end of for 'i' loop.

    } else {

	out <- "No matches found"
	attributes(out) <- a

    }

    class(out) <- "FeatureMatchAnalyzer"
    return(out)

} # end of 'FeatureMatchAnalyzer.matched.centmatch' function.

FeatureMatchAnalyzer.matched.deltamm <- function(x, which.comps=c("cent.dist", "angle.diff", "area.ratio", "int.area",
		    "bdelta", "haus", "ph", "med", "msd", "fom", "minsep", "bearing"),
    		    sizefac=1, alpha=0.1, k=4, p=2, c=Inf, distfun="distmapfun", ...,
    		    y=NULL, matches=NULL, object=NULL) {

    if(!is.null(matches)) obj <- matches
    else if(!is.null(y)) obj <- y
    else obj <- x

    loc <- attributes(obj)$loc

    Yfeats <- obj$Y.feats
    Xfeats <- obj$X.feats

    if(dim(obj$matches)[ 1 ] == 0) {
	# stop("FeatureAnalyzer: This function requires matches!")
	out <- "No matches found"
    } else {
        n <- dim(obj$matches)[1]

        out <- list()
        for(i in 1:n) out[[i]] <- FeatureComps(Y=Yfeats[[i]], X=Xfeats[[i]], which.comps=which.comps,
					sizefac=sizefac, alpha=alpha, k=k, p=p, c=c, distfun=distfun, loc = loc, ...)
    }

    class(out) <- "FeatureMatchAnalyzer"
    return(out)

} # end of 'FeatureMatchAnalyzer.matched.deltamm' function.

print.FeatureMatchAnalyzer <- function(x, ...) {

    if(is.list(x)) {

        n <- length(x)
        hold <- x[[1]]
        m <- length(hold)
        cn <- names(hold)

        for(i in 1:n) {

	    l <- unlist(lapply(x[[i]], length))
	    if(any(l==0)) x[[i]][[(1:m)[l==0]]] <- NA

        } # end of for 'i' loop.

        out <- c(unlist(x))
        attributes(out) <- NULL
        out <- matrix(out, n, m, byrow=TRUE)
        colnames(out) <- cn
        print(out)

        invisible(out)

    } else {

	attributes(x) <- NULL
	print(c(x))
	invisible()

    }

} # end of 'print.FeatureMatchAnalyzer' function.

summary.FeatureMatchAnalyzer <- function(object, ...) {

    if(is.list(object)) {

        args <- list(...)
        if(is.null(args$silent)) silent <- FALSE
        else silent <- args$silent
        n <- length(object)
        m <- length(object[[1]])
        res <- matrix(NA, n, m)
        colnames(res) <- names(object[[1]])

        for( i in 1:m) {

	    if(!silent) {

   	        if(names(object[[1]])[i] == "ph") cat("\n", "Partial Hausdorff Distance:\n")
   	        else if(names(object[[1]])[i] == "med") cat("\n", "Mean Error Distance:\n")
   	        else if(names(object[[1]])[i] == "msd") cat("\n", "Mean Square Error Distance:\n")
   	        else if(names(object[[1]])[i] == "fom") cat("\n", "Pratt\'s Figure of Merit:\n")
   	        else if(names(object[[1]])[i] == "minsep") cat("\n", "Minimum Separation Distance:\n")
   	        else if(names(object[[1]])[i] == "cent.dist") cat("\n", "Centroid Distance:\n")
   	        else if(names(object[[1]])[i] == "angle.diff") cat("\n", "Angle Difference:\n")
   	        else if(names(object[[1]])[i] == "area.ratio") cat("\n", "Area Ratio:\n")
   	        else if(names(object[[1]])[i] == "int.area") cat("\n", "Intersection Area:\n")
   	        else if(names(object[[1]])[i] == "bdelta") cat("\n", "Baddeley\'s Delta Metric:\n")
   	        else if(names(object[[1]])[i] == "haus") cat("\n", "Hausdorff Distance:\n")
   	        else if(names(object[[1]])[i] == "bearing") cat("\n", "Bearing:\n")

	    } # end of if '!silent' stmt.

	    hold <- numeric(n)+NA
	    for(j in 1:n) if(length(object[[j]][[i]])>0) hold[j] <- object[[j]][[i]]
	    res[,i] <- hold
	    if(!silent) print(hold)
       } # end of for 'i' loop.

   if(!is.null(args$interest)) {

	int <- args$interest
	a <- matrix(int, n, m, byrow=TRUE)

	if(!is.null(args$con)) {
	   con <- args$con
	   con <- match.fun(con)
	   a <- con(res, a, which.comps=names(object[[1]]))
	}

	out <- list()
	out$match.properties <- res
	b <- rowSums(a*res, na.rm=TRUE)
	out$object.interest <- b 
	out$interest.values <- int
	if(!silent) {
	   cat("\n", "Matched object interest values.\n")
	   print(b)
	}
   } else out <- res
       return(invisible(out))
    } else print(object)
    invisible()
} # end of 'summary.FeatureMatchAnalyzer' function.

plot.FeatureMatchAnalyzer <- function(x, ..., type=c("all", "ph", "med", "msd", "fom", "minsep",
					"cent.dist", "angle.diff", "area.ratio", "int.area", "bearing",
					"bdelta", "haus")) {

    if(is.list(x)) {

        type <- tolower(type)
        type <- match.arg(type)

        a <- attributes(x)
        y <- summary(x, silent=TRUE)
        n <- dim(y)[2]

	if(is.null(n)) {
	    n <- 1
	    y <- matrix(y, ncol=1)
	}

        if(type == "all") {

            if(!is.null(a$msg)) par(oma=c(0,0,2,0))
    
            for(i in 1:n) {

                if(is.element(colnames(y)[i], c("angle.diff","bearing"))) {

                    if(colnames(y)[i] == "angle.diff") t1 <- "Angle Difference"
                    else t1 <- "Bearing from Xhat to X\n obj. centroids (ref. = north)"
                    circ.plot(rad(y[,i]), main=t1, shrink=1.5)

    	        } else {

                    if(colnames(y)[i] == "ph") t1 <- "Partial Hausdorff \nDistance"
                    else if(colnames(y)[i] == "med") t1 <- "Mean Error Distance"
                    else if(colnames(y)[i] == "msd") t1 <- "Mean Square Error \nDistance"
                    else if(colnames(y)[i] == "fom") t1 <- "Pratt\'s Figure \nof Merit"
                    else if(colnames(y)[i] == "minsep") t1 <- "Minimum Separation \nDistance"
                    else if(colnames(y)[i] == "cent.dist") t1 <- "Centroid Distance"
                    else if(colnames(y)[i] == "area.ratio") t1 <- "Area Ratio"
                    else if(colnames(y)[i] == "int.area") t1 <- "Intersection Area"
                    else if(colnames(y)[i] == "bdelta") t1 <- "Baddeley\'s Delta Metric"
                    else if(colnames(y)[i] == "haus") t1 <- "Hausdorff Distance"
                    barplot(y[,i], main=t1, ...)

    	        }

            } # end of for 'i' loop.

        } else {

	    if(type=="ph") t1 <- "Partial Hausdorff \nDistance"
            else if(type == "med") t1 <- "Mean Error Distance"
            else if(type == "msd") t1 <- "Mean Square Error \nDistance"
            else if(type == "fom") t1 <- "Pratt\'s Figure \nof Merit"
            else if(type == "minsep") t1 <- "Minimum Separation \nDistance"
            else if(type == "cent.dist") t1 <- "Centroid Distance"
	    else if(type == "angle.diff") t1 <- "Angle Difference"
            else if(type == "area.ratio") t1 <- "Area Ratio"
            else if(type == "int.area") t1 <- "Intersection Area"
	    else if(type == "bearing") t1 <- "Bearing from Xhat to X\n obj. centroids (ref. = north)"
            else if(type == "bdelta") t1 <- "Baddeley\'s Delta Metric"
            else if(type == "haus") t1 <- "Hausdorff Distance"

	    i <- (1:ncol(y))[colnames(y) == type]

	    if(is.element(type, c("angle.diff", "bearing"))) circ.plot(rad(y[,i]), main=t1, shrink=1.5)
	    else barplot(y[,i], main=t1, ...)

        }

    } # end of if 'is.list(x)' stmts. 

    if(!is.null(a$msg)) mtext(a$msg, line=0.05, outer=TRUE)

    invisible()

} # end of 'plot.FeatureMatchAnalyzer' function.

bearing <- function(point1, point2, deg=TRUE, aty="compass") {

   if(is.null( dim( point1)) & length( point1)==2) point1 <- matrix(point1, 1, 2)
   if(is.null( dim( point2)) & length( point2)==2) point2 <- matrix(point2, 1, 2)
   
   if(deg) {

      # convert latitudes to radians
      point1[,2] <- point1[,2]*pi/180
      point2[,2] <- point2[,2]*pi/180

   } 

   # compute difference in longitude
   # dlon <- (point1[,1] - point2[,1])
   dlon <- (point1[,1] - point2[,1])
   if( deg) dlon <- dlon*pi/180
   S <- cos( point2[,2])*sin( dlon)
   Cval <- cos( point1[,2])*sin( point2[,2]) - sin( point1[,2])*cos( point2[,2])*cos( dlon)
   out <- atan2( S, Cval)
   
   # convert to degrees
   if( deg) out <- out*180/pi
   
   # out[out < 0 ] <-  out[out < 0] + 360
   # out[ out > 360] <- NA
   
   # print( out)
   
   if(aty == "radial") {

       # convert to polar coordinate type angle
       out[ (out >= 0)  &  (out <= 45)]  <- abs(out[ (out >= 0)  &  (out <= 45)] - 45) + 45
       out[ out > 45] <- 90 - out[ out > 45]
       out[ out < 0] <- out[ out < 0] + 360

      } # end of if aty is "radial" stmts.

   return( out)

} # end of 'bearing' function.

# compositer <- function(x, ...) {
#     UseMethod("compositer", x)
# } # end of 'compositer' function.
# 
# compositer.SpatialVx <- function(x, ..., time.point=NULL, model=1, identfun="threshsizer", verbose=FALSE) {
# 
#     if(verbose) begin.tiid <- Sys.time()
# 
#     theCall <- match.call()
# 
#     a <- attributes(x)
# 
#     ## Begin: Get the data sets
# #     if(!is.null(time.point)) {
# #         if(!missing(time.point) && !missing(model)) dat <- datagrabber(object, time.point=time.point, model=model)
# #         else if(!missing(time.point)) dat <- datagrabber(object, time.point=time.point)
# #         else if(!missing(model)) dat <- datagrabber(object, model=model)
# #         else dat <- datagrabber(object)
# # 
# #         X <- dat$X
# #         Xhat <- dat$Xhat
# #     } else {
# # 	X <- x[[1]]
# # 	if(is.list(x[[2]])) {
# # 	    if(is.function(model)) Xhat <- do.call(model, c(list(x[[2]]), list(...)))
# # 	    else if(length(model) == 1) Xhat <- x[[2]][[model]]
# # 	    else stop("compositer.SpatialVx: invalid model argument.")
# # 	} else Xhat <- x[[2]]
# #     }
#     ## End: Get the data sets
# 
#     ## Internal function to apply compositing strategy to each
#     ## forecast model and to verification field(s).
#     cfun <- function(Z, ..., ifun, time.point, model) {
# 
# 	o <- do.call(ifun, c(list(object=Z, time.point=time.point, model=model), list(...)))
# 
# 	ff <- function(x, obj) {
# 	    y <- c(obj)
# 	    z <- as.logical(c(x$m))
# 	    y[ !z ] <- 0
# 	    return(y)
# 	} # end of internal-internal 'ff' function.
# 
# 	X <- lapply(o$X.feats, ff, obj=o$X)
# 	Xhat <- lapply(o$Y.feats, ff, obj=o$Xhat)
# 
# 	nobj1 <- length(X)
# 	nobj2 <- length(Xhat)
# 
# 	out <- list()
# 	out$X <- matrix(unlist(X), nrow=nobj1, byrow=TRUE)
# 	out$Xhat <- matrix(unlist(Xhat), nrow=nobj2, byrow=TRUE)
# 
# 	return(out)
#     } # end of internal 'cfun' function.
# 
#     out <- list()
#     attributes(out) <- a
# 
#     loop.through.time <- (length(a$xdim) > 2) && (is.null(time.point) || length(time.point) > 1)
#     if(loop.through.time) {
# 	if(is.character(time.point)) tlab <- time.point
#         else tlab <- NULL
#     }    
# 
#     if(is.null(time.point) && length(a$xdim) == 3) time.point <- 1:a$xdim[3]
#     else if(is.null(time.point) && length(a$xdim) == 2) time.point <- 1
#     if(!is.numeric(time.point)) time.point <- (1:a$xdim[3])[time.point == a$time]
# 
#     if(!loop.through.time) res <- cfun(Z=x, ..., time.point=time.point, model=model, ifun=identfun)
#     else {
# 
# 	nvx <- dim(X)[3]
# 	res <- list()
# 	res$X <- numeric(0)
# 	res$Xhat <- numeric(0)
# 
# 	if(verbose) cat("\n", "Looping through time.\n")
# 	for(i in 1:nvx) {
# 	    if(verbose) cat(time.point[i], "\n")
# 	    tmp <- cfun(x, ..., fun=identfun, time.point=time.point[i], model=model)
# 	    X <- tmp$X
# 	    Xhat <- tmp$Xhat
# 
# 	    n1 <- dim(X)
# 	    if(is.null(n1)) n1 <- 1
# 
# 	    n2 <- dim(Xhat)
# 	    if(is.null(n2)) n2 <- 1
# 
# 	    if(!is.null(X)) rownames(X) <- paste(tlab[i], ".obj", 1:n1, sep="") 
# 	    if(!is.null(Xhat)) rownames(Xhat) <- paste(tlab[i], ".obj", 1:n2, sep="")
# 
# 	    res$X <- rbind(res$X, X)
# 	    res$Xhat <- rbind(res$Xhat, Xhat)
# 
# 	} # end of for 'i' loop.
# 
#     } # end of if else more than one time point stmts.
# 
#     # TO DO: Should now have matrices with rows corresponding to unique objects in each
#     # field.  Need to check this.  Next, need to center each object onto a new relative
#     # grid.  Will not worry about the relative grid size here.  Will allow for that in
#     # subsequent functions.
#     # out$X.composites <- res$X
#     # out$Xhat.composites <- res$Xhat
# 
#     attr(out, "time.point") <- time.point
#     attr(out, "model") <- model
#     attr(out, "call") <- theCall
# 
#     attr(out, "identifier.function") <- identfun
# 
#     if(verbose) print(Sys.time() - begin.tiid)
# 
#     class(out) <- "compositer"
#     return(out)
# } # end of 'compositer.SpatialVx' function.
# 
# compositer.default <- function(x, ..., xhat, identfun="threshsizer") {
#     tmp <- make.SpatialVx(x, xhat, data.name=c(deparse(substitute(x)), deparse(substitute(xhat))))
#     res <- compositer.SpatialVx(tmp, ..., identfun=identfun)
#     return(res)
# } # end of 'compositer.default' function.
