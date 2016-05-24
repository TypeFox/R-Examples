compositer <- function(x, level = 0, verbose = FALSE, ...) {

    if(!is.element(class(x), c("features", "matched", "combined"))) stop("compositer: invalid x argument.")
    UseMethod("compositer", x)

} # end of 'compositer' function.

compositer.matched <- function(x, level = 0, verbose = FALSE, ...) {

    class(x) <- "features"
    UseMethod("compositer", x)

} # end of 'compositer.matched' function.

compositer.combined <- function(x, level = 0, verbose = FALSE, ...) {

    if(class(x) != "combined") stop("compositer.combined: invalid x argument.")

    if(verbose) begin.tiid <- Sys.time()

    if(is.null(x$X.feats) || length(x$X.feats) == 0 || is.null(x$Y.feats) || length(x$Y.feats) == 0) {

        if(verbose) cat("No features in or both fields.  Returning NULL.")
        return(NULL)

    }

    out <- x

    xdim <- attributes(x)$xdim
    if(is.null(xdim)) {

	if(length(xdim) == 3) xdim <- dim(x$X[,,1])
	else if(length(xdim) == 2) xdim <- dim(x$X)
	else stop("compositer.combined: incorrect dimensions for X.")

    }

    loc <- cbind(rep(1:xdim[1], xdim[2]), rep(1:xdim[2], each = xdim[1]))

    # Approximating the centroid to nearest grid square.  Might be better to interpolate the
    # field(s), but at greater computational cost.
    cen <- round(colMeans(loc), digits = 0)

    ifun <- function(y, loc) {

        res <- FeatureProps(x = y, which.props = "centroid", loc = loc)
        res <- c(res$centroid$x, res$centroid$y)
        names(res) <- c("x", "y")

        return(res)

    } # end of internal 'ifun' internal function.

    Xcenters <- round(matrix(unlist(lapply(x$X.feats, ifun, loc = loc)), ncol = 2, byrow = TRUE), digits = 0)
    Ycenters <- round(matrix(unlist(lapply(x$Y.feats, ifun, loc = loc)), ncol = 2, byrow = TRUE), digits = 0)

    # calcualte the centroid distances between each pair of features.
    A <- outer(Xcenters[,1], Ycenters[,1], "-")
    B <- outer(Xcenters[,2], Ycenters[,2], "-")
    d <- sqrt(A^2 + B^2)

    # find the minimum centroid distance between each feature in one field and all features in the other.
    # 'dX' is min. centroid distance between each obs feature and all forecast features.
    # 'dY' is min. centroid distance between each forecast feature and all obs features.
    dX <- apply(d, 1, min, na.rm = TRUE)
    dY <- apply(d, 2, min, na.rm = TRUE)

    out$distances <- list(X = dX, Xhat = dY)

    Xtrans <- matrix(cen, nrow = nrow(Xcenters), ncol = 2, byrow = TRUE) - Xcenters
    Ytrans <- matrix(cen, nrow = nrow(Ycenters), ncol = 2, byrow = TRUE) - Ycenters

    # Find new dimensions and locations such that they range from 1 to K,
    # and each feature will be included in its entirety in the new locations.
    # These new coordinates may be too large, but we will tighten them up later.

    newdim1 <- xdim[1] + max(abs(c(Xtrans[,1], Ytrans[,1])), na.rm = TRUE)
    newdim2 <- xdim[2] + max(abs(c(Xtrans[,2], Ytrans[,2])), na.rm = TRUE)
    newdim <- c(newdim1, newdim2)

    newloc <- cbind(rep(1:newdim[1], newdim[2]), rep(1:newdim[2], each = newdim[1]))
    newcen <- round(colMeans(newloc, na.rm = TRUE), digits = 0)

    # Now, figure out the translation required for the new grid.
    newtr <- newcen - cen
    Xtrans <- Xtrans + matrix(newtr, nrow = nrow(Xtrans), ncol = 2, byrow = TRUE)
    Ytrans <- Ytrans + matrix(newtr, nrow = nrow(Ytrans), ncol = 2, byrow = TRUE)

    mat <- matrix(0, newdim[1], newdim[2])

    tfun <- function(id, obj, loc0, loc1, trans, mat, ff) {

        Im <- obj[[ id[1] ]]
        ind <- c(Im$m)

        if(!any(ind)) return(NULL)

        ind0 <- loc0[ind, , drop = FALSE]
        ind1 <- ind0 + matrix(trans[id[1], ], nrow = nrow(ind0), ncol = 2, byrow = TRUE)
        # mat[ ind1 ] <- 1
	ff2 <- ff[,, id[2]]
        mat[ ind1 ] <- ff2[ Im$m ]

	return(mat)

    } # end of internal 'tfun' function.

    indX <- as.list(1:dim(Xtrans)[1])
    indY <- as.list(1:dim(Ytrans)[1])

    indX2 <- apply(x$X.labeled, 3, max, na.rm = TRUE)
    indY2 <- apply(x$Y.labeled, 3, max, na.rm = TRUE)

    indX3 <- indY3 <- numeric(0)
    for(i in 1:length(indX2)) indX3 <- c(indX3, rep(i, indX2[ i ]))
    for(i in 1:length(indY2)) indY3 <- c(indY3, rep(i, indY2[ i ]))

    for(i in 1:length(indX)) indX[[ i ]] <- c(indX[[ i ]], indX3[ i ])
    for(i in 1:length(indY)) indY[[ i ]] <- c(indY[[ i ]], indY3[ i ])

    # Results below should be a new list of lists of big matrices containing each feature as before,
    # but where each has the same centroid and is on a larger grid.
    resX <- lapply(indX, tfun, obj = x$X.feats, loc0 = loc, loc1 = newloc, trans = Xtrans, mat = mat, ff = x$X)
    resY <- lapply(indY, tfun, obj = x$Y.feats, loc0 = loc, loc1 = newloc, trans = Ytrans, mat = mat, ff = x$Xhat)

    # Now try to find the smallest grid that contains each feature, and return as a component of 'out'.

    rg <- shrinkgrid(obj1 = resX, obj2 = resY, loc = newloc, level = level, verbose = verbose)

    out$Xcentered <- rg$obj1
    out$Ycentered <- rg$obj2

    if(verbose) print(Sys.time() - begin.tiid)
    class(out) <- "composited"
    return(out)

} # end of 'compositer.combined' function.

compositer.features <- function(x, level = 0, verbose = FALSE, ...) {

    if(!is.element(class(x), c("features", "matched"))) stop("compositer: invalid x argument.")

    if(verbose) begin.tiid <- Sys.time()

    if(is.null(x$X.feats) || length(x$X.feats) == 0 || is.null(x$Y.feats) || length(x$Y.feats) == 0) {

	if(verbose) cat("No features in or both fields.  Returning NULL.")
	return(NULL)

    }

    out <- x

    xdim <- attributes(x)$xdim
    if(is.null(xdim)) xdim <- dim(x$X)

    loc <- cbind(rep(1:xdim[1], xdim[2]), rep(1:xdim[2], each = xdim[1]))

    # Approximating the centroid to nearest grid square.  Might be better to interpolate the
    # field(s), but at greater computational cost.
    cen <- round(colMeans(loc), digits = 0)

    ifun <- function(y, loc) {

	res <- FeatureProps(x = y, which.props = "centroid", loc = loc)
	res <- c(res$centroid$x, res$centroid$y)
	names(res) <- c("x", "y")

	return(res)

    } # end of internal 'ifun' internal function.

    Xcenters <- round(matrix(unlist(lapply(x$X.feats, ifun, loc = loc)), ncol = 2, byrow = TRUE), digits = 0)
    Ycenters <- round(matrix(unlist(lapply(x$Y.feats, ifun, loc = loc)), ncol = 2, byrow = TRUE), digits = 0)

    # calcualte the centroid distances between each pair of features.
    A <- outer(Xcenters[,1], Ycenters[,1], "-")
    B <- outer(Xcenters[,2], Ycenters[,2], "-")
    d <- sqrt(A^2 + B^2)

    # find the minimum centroid distance between each feature in one field and all features in the other.
    # 'dX' is min. centroid distance between each obs feature and all forecast features.
    # 'dY' is min. centroid distance between each forecast feature and all obs features.
    dX <- apply(d, 1, min, na.rm = TRUE)
    dY <- apply(d, 2, min, na.rm = TRUE)

    out$distances <- list(X = dX, Xhat = dY)

    Xtrans <- matrix(cen, nrow = nrow(Xcenters), ncol = 2, byrow = TRUE) - Xcenters
    Ytrans <- matrix(cen, nrow = nrow(Ycenters), ncol = 2, byrow = TRUE) - Ycenters


    # Find new dimensions and locations such that they range from 1 to K,
    # and each feature will be included in its entirety in the new locations.
    # These new coordinates may be too large, but we will tighten them up later.

    newdim1 <- xdim[1] + max(abs(c(Xtrans[,1], Ytrans[,1])), na.rm = TRUE)
    newdim2 <- xdim[2] + max(abs(c(Xtrans[,2], Ytrans[,2])), na.rm = TRUE) 
    newdim <- c(newdim1, newdim2)

    newloc <- cbind(rep(1:newdim[1], newdim[2]), rep(1:newdim[2], each = newdim[1]))
    newcen <- round(colMeans(newloc, na.rm = TRUE), digits = 0)

    # Now, figure out the translation required for the new grid.
    newtr <- newcen - cen
    Xtrans <- Xtrans + matrix(newtr, nrow = nrow(Xtrans), ncol = 2, byrow = TRUE)
    Ytrans <- Ytrans + matrix(newtr, nrow = nrow(Ytrans), ncol = 2, byrow = TRUE)

    mat <- matrix(0, newdim[1], newdim[2])

    tfun <- function(id, obj, loc0, loc1, trans, mat, ff) {

	Im <- obj[[ id ]]
	ind <- c(Im$m)

	if(!any(ind)) return(NULL)

	ind0 <- loc0[ind, , drop = FALSE]
	ind1 <- ind0 + matrix(trans[id, ], nrow = nrow(ind0), ncol = 2, byrow = TRUE)
	# mat[ ind1 ] <- 1
	mat[ ind1 ] <- ff[ Im$m ]

	return(mat)

    } # end of internal 'tfun' function.

    indX <- as.list(1:dim(Xtrans)[1])
    indY <- as.list(1:dim(Ytrans)[1])

    # Results below should be a new list of lists of owin objects containing each feature as before,
    # but where each has the same centroid and is on a larger grid.
    resX <- lapply(indX, tfun, obj = x$X.feats, loc0 = loc, loc1 = newloc, trans = Xtrans, mat = mat, ff = x$X)
    resY <- lapply(indY, tfun, obj = x$Y.feats, loc0 = loc, loc1 = newloc, trans = Ytrans, mat = mat, ff = x$Xhat)

    # Now try to find the smallest grid that contains each feature, and return as a component of 'out'.

    rg <- shrinkgrid(obj1 = resX, obj2 = resY, loc = newloc, level = level, verbose = verbose)

    out$Xcentered <- rg$obj1
    out$Ycentered <- rg$obj2

    if(verbose) print(Sys.time() - begin.tiid)

    class(out) <- "composited"
    return(out)

} # end of 'compositer.features' function.

shrinkgrid <- function(obj1, obj2, loc, level, verbose) {

    # Try to find the smallest grid that contains each feature.

    bfun <- function(x, loc, lvl, vb) {

        id <- x > lvl
        l <- loc[ c(id), , drop = FALSE]
        res <- c(range(l[,1], finite = TRUE), range(l[,2], finite = TRUE))
        names(res) <- c("x1", "x2", "y1", "y2")
        if(vb) print(res)
        return(res)

    } # end of internal 'bfun' function.

    if(verbose) cat("\n", "Finding the smallest box to make the relative grid.\n")
    crnsX <- lapply(obj1, bfun, loc = loc, lvl = level, vb = verbose)
    crnsY <- lapply(obj2, bfun, loc = loc, lvl = level, vb = verbose)

    grabby <- function(x, l) return(x[l])

    Xx1 <- min(c(unlist(lapply(crnsX, grabby, l = "x1"))), na.rm = TRUE)
    Xx2 <- max(c(unlist(lapply(crnsX, grabby, l = "x2"))), na.rm = TRUE)
    Xy1 <- min(c(unlist(lapply(crnsX, grabby, l = "y1"))), na.rm = TRUE)
    Xy2 <- max(c(unlist(lapply(crnsX, grabby, l = "y2"))), na.rm = TRUE)

    Yx1 <- min(c(unlist(lapply(crnsY, grabby, l = "x1"))), na.rm = TRUE)
    Yx2 <- max(c(unlist(lapply(crnsY, grabby, l = "x2"))), na.rm = TRUE)
    Yy1 <- min(c(unlist(lapply(crnsY, grabby, l = "y1"))), na.rm = TRUE)
    Yy2 <- max(c(unlist(lapply(crnsY, grabby, l = "y2"))), na.rm = TRUE)

    rc <- cbind(c(min(Xx1, Yx1, na.rm = TRUE), max(Xx2, Yx2, na.rm = TRUE)),
		c(min(Xy1, Yy1, na.rm = TRUE), max(Xy2, Yy2, na.rm = TRUE)))

    if(verbose) {

        rgc <- rc
        colnames(rgc) <- c("x", "y")
        cat("Relative grid corners:\n")
        print(rgc)

    }

    rebfun <- function(x, cnr) {

        xd <- dim(x)

        if(cnr[1,1] < 1) cnr[1,1] <- 1
        if(cnr[2,1] > xd[1]) cnr[2,1] <- xd[1]
        if(cnr[1,2] < 1) cnr[1,2] <- 1
        if(cnr[2,2] > xd[2]) cnr[2,2] <- xd[2]

        out <- x[ cnr[1,1]:cnr[2,1], cnr[1,2]:cnr[2,2] ]
        return(out)

    } # end of internal 'rebfun' function.

    obj1 <- lapply(obj1, rebfun, cnr = rc)
    obj2 <- lapply(obj2, rebfun, cnr = rc)

    # obj1 <- lapply(obj1, im)
    # obj2 <- lapply(obj2, im)

    return(list(obj1 = obj1, obj2 = obj2))

} # end of 'shrinkgrid' function.

plot.composited <- function(x, ..., type = c("all", "X", "Xhat", "X|Xhat", "Xhat|X"), dist.crit = 100, FUN = "mean",
    col = c("gray", tim.colors(64)) ) {

    type <- match.arg(type)

    if(is.element(type, c("all", "X|Xhat"))) {

	dX <- x$distances$X
	if(all(dX > dist.crit)) {

	    if(type == "X|Xhat") {

		warning(paste("plot: No minimum centroid distances less than distance criterion = ", dist.crit, ".  Nothing to plot.", sep = ""))
		invisible()

	    } else warning(paste("plot: No minimum centroid distances less than distance criterion = ", dist.crit, 
			".  Nothing to plot for X|Xhat.", sep = ""))
	}

    }

    if(is.element(type, c("all", "Xhat|X"))) {

	dY <- x$distances$Xhat

	if(all(dY > dist.crit)) {

            if(type == "Xhat|X") {

                warning(paste("plot: No minimum centroid distances less than distance criterion = ", dist.crit, ".  Nothing to plot.", sep = ""))
                invisible()

            } else warning(paste("plot: No minimum centroid distances less than distance criterion = ", dist.crit, 
			".  Nothing to plot for Xhat|X.", sep = ""))
        }

    }

    a <- attributes(x)

    if(!is.null(a$data.name)) {

        dn <- a$data.name
        if(length(dn) == 3) {

	    msg <- dn[1]
	    dn <- dn[-1]
	    
        } else msg <- NULL

	Xname <- dn[1]
        Xhat.name <- dn[2]
        X.Xhat.name <- paste(dn[1], dn[2], sep = "|")
        Xhat.X.name <- paste(dn[2], dn[1], sep = "|")

    } else {

        Xname <- "X"
        Xhat.name <- "Xhat"
        X.Xhat.name <- "X|Xhat"
        Xhat.X.name <- "Xhat|X"
	msg <- NULL

    }

    if(!is.null(msg)) par(oma = c(0, 0, 2, 0))

    n <- length(x$Xcentered)
    m <- length(x$Ycentered)

    if(m > 0) {

        if(is.element(type, c("all", "Xhat", "Xhat|X"))) {

	    # xd <- dim(as.matrix.im(x$Ycentered[[ 1 ]]))
	    xd <- dim(x$Ycentered[[ 1 ]])

	    if(is.element(type, c("all", "Xhat"))) {

	        # Xhat <- as.matrix.im(x$Ycentered[[ 1 ]])
		Xhat <- x$Ycentered[[ 1 ]]
	        # if(m > 1) for(i in 2:m) Xhat <- array(c(c(Xhat), c(as.matrix.im(x$Ycentered[[ i ]]))), dim = c(xd, i))
		 if(m > 1) for(i in 2:m) Xhat <- array(c(c(Xhat), c(x$Ycentered[[ i ]])), dim = c(xd, i))
	        # Xhat <- t(apply(Xhat, 1:2, FUN))
		Xhat <- apply(Xhat, 1:2, FUN)
	        Z <- Xhat
	        zl1 <- range(c(Xhat), finite = TRUE)

	    } else zl1 <- Inf

	    if(is.element(type, c("all", "Xhat|X"))) {

	        if(any(dY <= dist.crit)) {

	            XhatX0 <- x$Ycentered[ dY <= dist.crit ]
	            # XhatX <- as.matrix.im(XhatX0[[ 1 ]])
		    XhatX <- XhatX0[[ 1 ]]
		    # if(length(XhatX0) > 1) for(i in 2:length(XhatX0)) XhatX <- array(c(c(XhatX), c(as.matrix.im(XhatX0[[ i ]]))), dim = c(xd, i) )
		    if(length(XhatX0) > 1) for(i in 2:length(XhatX0)) XhatX <- array(c(c(XhatX), c(XhatX0[[ i ]])), dim = c(xd, i))
		    # XhatX <- t(apply(XhatX, 1:2, FUN))
		    XhatX <- apply(XhatX, 1:2, FUN)
		    Z <- XhatX

	        } else XhatX <- x$Ycentered[[ 1 ]] * 0
# XhatX <- t(as.matrix.im(x$Ycentered[[ 1 ]])) * 0

	        zl1 <- range(c(zl1, c(XhatX)), finite = TRUE)

	    }

        } else zl1 <- Inf

    } else if(is.element(type, c("all", "Xhat", "Xhat|X"))) {

	Xhat <- Xhat0 <- matrix(0, a$xdim[1], a$xdim[2])
	zl1 <- c(0, 0)

    }

    if(n > 0) {

        if(is.element(type, c("all", "X", "X|Xhat"))) {

	    # xd2 <- dim(as.matrix.im(x$Xcentered[[ 1 ]]))
	    xd2 <- dim(x$Xcentered[[ 1 ]])

	    if(is.element(type, c("all", "X"))) {

	        # X <- as.matrix.im(x$Xcentered[[ 1 ]])
		X <- x$Xcentered[[ 1 ]]
	        # if(n > 1) for(i in 2:n) X <- array(c(c(X), c(as.matrix.im(x$Xcentered[[ i ]]))), dim = c(xd2, i))
		# X <- t(apply(X, 1:2, FUN))
		if(n > 1) for(i in 2:n) X <- array(c(c(X), c(x$Xcentered[[ i ]])), dim = c(xd2, i))
		X <- apply(X, 1:2, FUN)
		Z <- X
		zl2 <- range(c(X), finite = TRUE)

	    } else zl2 <- Inf

	    if(is.element(type, c("all", "X|Xhat"))) {

		if(any(dX <= dist.crit)) {

		    XXhat0 <- x$Xcentered[ dX <= dist.crit ]
		    # XXhat <- as.matrix.im(XXhat0[[ 1 ]])
		    XXhat <- XXhat0[[ 1 ]]
		    # if(length(XXhat0) > 1) for(i in 2:length(XXhat0)) XXhat <- array(c(c(XXhat), c(as.matrix.im(x$Xcentered[[ i ]]))), dim = c(xd2, i))
		    # XXhat <- t(apply(XXhat, 1:2, FUN))
		    if(length(XXhat0) > 1) for(i in 2:length(XXhat0)) XXhat <- array(c(c(XXhat), c(x$Xcentered[[ i ]])), dim = c(xd2, i))
		    XXhat <- apply(XXhat, 1:2, FUN)
		    Z <- XXhat
		    zl2 <- range(c(zl2, c(XXhat)), finite = TRUE)

		} else XXhat <- x$Xcentered[[ 1 ]] * 0
# XXhat <- t(as.matrix.im(x$Xcentered[[ 1 ]])) * 0

	    }

        }

    } else {

	X <- XXhat <- matrix(0, a$xdim[1], a$xdim[2])
	zl2 <- c(0, 0)

    }

    if(type == "all") {

	par(mfrow = c(2, 2))
	zl <- range(c(zl1, zl2), finite = TRUE)

    } else if(type == "X|Xhat") zl <- zl2
    else if(type == "Xhat|X") zl <- zl1

    if(is.element(type, c("all", "X"))) image(X, col = col, zlim = zl, main = Xname)
    if(is.element(type, c("all", "Xhat"))) image(Xhat, col = col, zlim = zl, main = Xhat.name)
    if(is.element(type, c("all", "X|Xhat"))) image(XXhat, col = col, zlim = zl, main = X.Xhat.name)
    if(is.element(type, c("all", "Xhat|X"))) image(XhatX, col = col, zlim = zl, main = Xhat.X.name)
    image.plot(Z, col = col, zlim = zl, legend.only = TRUE, ...)

    if(!is.null(msg)) {

	title("")
	mtext(paste(msg, " (composite ", FUN, ")", sep = ""), line = 0.05, outer = TRUE)

    }

    invisible()

} # end of 'plot.composited' function.
