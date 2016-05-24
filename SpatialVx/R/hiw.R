hiw <- function(x, simplify = 0, A = pi * c(0, 1/16, 1/8, 1/6, 1/4, 1/2, 9/16, 5/8, 2/3, 3/4), verbose = FALSE, ...) {

    if(any(A < 0) || any(A > pi)) stop("hiw: A must give angles in radians between 0 and pi.")

    dupang <- (A + pi) == A
    if(any(dupang)) A <- A[ !dupang ]

    obs <- x$X.feats
    mod <- x$Y.feats

    xdim <- dim(x$X)

    if(simplify < 0 || is.infinite(simplify)) stop("hiw: invalid simplify argument.")

    out <- x

    if(verbose) cat("Converting features to polygon objects.\n")
    
    if(is.na(simplify) || simplify == 0) {

	opoly <- lapply(obs, as.polygonal)
	mpoly <- lapply(mod, as.polygonal)

    } else {

	pfun <- function(x,s) return(as.polygonal(simplify.owin(x, s)))
	opoly <- lapply(obs, pfun, s = simplify)
	mpoly <- lapply(mod, pfun, s = simplify)

    } # end of if else 'simplify' object before doing shape analysis stmts.

    if(verbose) cat("Polygon objects created.\n")

    # opoly <- lapply(opoly, as.psp)
    # mpoly <- lapply(mpoly, as.psp)
    if(verbose) cat("Finding the edges of each feature in each field.\n")
    opoly <- lapply(opoly, edges)
    mpoly <- lapply(mpoly, edges)
    if(verbose) cat("Edges found.\n")

    if(verbose) cat("finding the centroids of each feature in each field.\n")
    ocen <- lapply(obs, centroid.owin)
    mcen <- lapply(mod, centroid.owin)
    if(verbose) cat("Centroids found.\n")

    ocen <- cbind(c(unlist(lapply(ocen, function(x) return(x$x)))),
		c(unlist(lapply(ocen, function(x) return(x$y)))))

    mcen <- cbind(c(unlist(lapply(mcen, function(x) return(x$x)))),
		c(unlist(lapply(mcen, function(x) return(x$y)))))

    d <- prod(xdim)
    ol <- rep(d, length(obs))
    ml <- rep(d, length(mod))

    ifun <- function(id, x, xmid, ymid, l, a, vb) {

	if(vb) cat(id, " ")

	X <- x[[ id ]]
	N <- length(a)

	# At most, there should be twice the number of input angles of feature-crossing points
        # per feature.
        numseg <- 2 * N

	# Set up a line segment object with the user chosen angles such that the lengths stretch beyond the
	# extent of the feature in order to ensure they will overlap with the feature (may be a problem if the feature
	# is actually a cluster of features, but not going to worry about it).

	seg <- as.psp(list(xmid = rep(xmid[ id ], N), ymid = rep(ymid[ id ], N), length = rep(l[ id ], N), angle = a),
		window = as.rectangle(X))

	# find where the lines cross the boundary of the features (only the feature boundaries have been passed in).

	crs <- crossing.psp(seg, X)

	# Re-create the line segment objects so that the lines only stretch from the centroid to the boundaries.
	# May have multiple crossings for non-convex features and/or features with holes inside them.

	res <- as.psp(list(x0 = rep(xmid[ id ], crs$n), y0 = rep(ymid[ id ], crs$n),
			x1 = crs$x, y1 = crs$y), window = as.rectangle(X))

	# Try to remove unwanted line segments (e.g., if there is a hole, there could be several
	# unwanted line segments along the same angle).  Also, figure out if any line segments:
	# (i) cross only on one side of the mid-point, (ii) cross in one or fewer places.
	a2 <- round(angles.psp(res), digits = 6)
	l2 <- lengths.psp(res)

	mat <- res$ends
	dx <- zapsmall(mat$x1 - mat$x0) # + right, - left, 0 vertical
	dy <- zapsmall(mat$y1 - mat$y0) # + up, - down, 0 horizontal

	o <- order(a2)
        a2 <- a2[ o ]
	l2 <- l2[ o ]
	dx <- dx[ o ]
	dy <- dy[ o ]
        mat2 <- mat[o, ]

	m <- res$n

	QI   <- dx > 0 & dy > 0
	QII  <- dx < 0 & dy > 0
	QIII <- dx < 0 & dy < 0
	QIV  <- dx > 0 & dy < 0
	vert <- dx == 0 & dy != 0
	horz <- dx != 0 & dy == 0
	sing <- dx == 0 & dy == 0
	if(any(sing)) warning(paste("single point at angle: ", a2, " removing/ignoring.", sep = ""))

	# for each unique angle, need to check where the points lie.
	ua <- unique(a2)
	nua <- length(ua)

	outvec <- numeric(numseg) + NA

	sider <- function(qi, qii, qiii, qiv, v, h, Dx, Dy, len) {

	    # Once there are only two crossings, return a length 2 vector 
	    # with 1 or -1 depending on whether they are on opposite sides
	    # of each other (both 1), or the same side (shorter lengthed
	    # segment will be -1 and longer lengthed segment will be 1.

	    out <- numeric(2) + NA

	    if((any(qi) & any(qiii)) || (any(qii) & any(qiv))) {

                # Boundary crossings are on opposite sides.  Nothing 
                # more to do, except return out with 1's.
                out <- c(1, 1)

            } else if(all(qi) || all(qii) || all(qiii) || all(qiv)) {

                # Boundary crossings are all on the same side.
		# return a -1 for the closer point and +1 for the farther point.

                ind2 <- len == min(len) 
		out <- numeric(2)
                out[  ind2 ] <- -1
                out[ !ind2 ] <- 1

            } else if(all(v)) {

                # Vertical line situation.  
                Dy <- dy[ ind ]
		if(length(Dy) == 1) out[1] <- 1
		else if((Dy[1] > 0 & Dy[2] < 0) || (Dy[1] < 0 & Dy[2] > 0)) out <- c(1, 1)
                else {

                    ind2 <- len == min(len)
                    out[  ind2 ] <- -1
                    out[ !ind2 ] <- 1

                }

            } else if(all(h)) {

                Dx <- dx[ ind ]
		if(length(Dx) == 1) out[1] <- 1
                else if((Dx[1] > 0 & Dx[2] < 0) || (Dx[1] < 0 & Dx[2] > 0)) out <- c(1, 1)
                else {

                    ind2 <- len == min(len)
                    out[  ind2 ] <- -1
                    out[ !ind2 ] <- 1

                }

            } # end of if else opposite or same side stmts.

	    return(out)

	} # end of internal internal 'sider' function.

	extra.points <- integer(0)

	# These are the angles that should be there.
        a3 <- rep(round(a, digits = 6), each = 2)

	for(i in 1:nua) {

	    ind <- (1:m)[ a2 == ua[ i ] ]
	    bigind <- (1:numseg)[ a3 == ua[ i ] ]
	    num <- length(ind)

            sg <- sing[ ind ]

	    if(any(sg)) warning("hiw: single point.  Don't know what to do with it.")

	    if(num == 2) {

		# Simplest case.  Need only determine whether the
                # crossings are on the same side or opposite sides
                # of the centroid.  If opposite, indicate with a +1
                # in the returned shape attribute vector (otherwise
                # indicate with a -1).

		qi <- QI[ ind ]
                qii <- QII[ ind ]
                qiii <- QIII[ ind ]
                qiv <- QIV[ ind ]
                v <- vert[ ind ]
                h <- horz[ ind ]

		outvec[ bigind ] <- sider(qi = qi, qii = qii, qiii = qiii, qiv = qiv, v = v, h = h, Dx = dx[ ind ], Dy = dy[ ind ], len = l2[ ind ])

	    } else if(num > 2) {

		# Must remove extra crossings segments before repeating above exercise.

		l3 <- l2[ ind ]
		qi <- QI[ ind ]
            	qii <- QII[ ind ]
            	qiii <- QIII[ ind ]
            	qiv <- QIV[ ind ]
            	v <- vert[ ind ]
            	h <- horz[ ind ]

		if(!any((any(qi) & any(qiii)) | (any(qii) & any(qiv)) | v | h)) {

		    # crossings are not on opposite sides of the centroid.

		    lmin <- l3 == min(l3, na.rm = TRUE)
                    lmax <- l3 == max(l3, na.rm = TRUE)

		    ind2 <- ind[ lmin ][ 1 ]
		    ind2 <- c(ind2, ind[ lmax ][ 1 ])

		} else {

		    # crossings are on opposite sides of the centroid.

		    if(any(qi) & any(qiii)) {

			lmax1 <- max(l3[ qi ], na.rm = TRUE)
			lmax2 <- max(l3[ qiii ], na.rm = TRUE)

			ind2 <- c(ind[ (l3 == lmax1) & qi ][ 1 ], ind[ (l3 == lmax2) & qiii ][ 1 ])

		    } else if(any(qii) & any(qiv)) {

			lmax1 <- max(l3[ qii ], na.rm = TRUE)
			lmax2 <- max(l3[ qiv ], na.rm = TRUE)

			ind2 <- c(ind[ (l3 == lmax1) & qii ][ 1 ], ind[ (l3 == lmax2) & qiv ][ 1 ])

		    } else if(any(v)) {

			upside <- dy[ ind ] > 0
			downside <- dy[ ind ] < 0
			lmax1 <- max(l3[ upside ], na.rm = TRUE)
			lmax2 <- max(l3[ downside ], na.rm = TRUE)

			ind2 <- c(ind[ (l3 == lmax1) & upside ][ 1 ], ind[ (l3 == lmax2) & downside ][ 1 ])

		    } else if(any(h)) {

			rightside <- dx[ ind ] > 0
			leftside <- dx[ ind ] < 0
			lmax1 <- max(l3[ rightside ], na.rm = TRUE)
			lmax2 <- max(l3[ leftside ], na.rm = TRUE)

			ind2 <- c(ind[ (l3 == lmax1) & rightside ][ 1 ], ind[ (l3 == lmax2) & leftside ][ 1 ])

		    }

		} # end of if crossings are on same or opposite sides stmts.

		extra.points <- c(extra.points, ind[ !is.element(ind, ind2) ])
		ind <- ind[ is.element(ind, ind2) ]

                qi <- QI[ ind ]
                qii <- QII[ ind ]
                qiii <- QIII[ ind ]
                qiv <- QIV[ ind ]
                v <- vert[ ind ]
                h <- horz[ ind ]

		outvec[ bigind ] <- sider(qi = qi, qii = qii, qiii = qiii, qiv = qiv, v = v, h = h, Dx = dx[ ind ], Dy = dy[ ind ], len = l2[ ind ])


	    } else if(num == 1) {

		outvec[ bigind ] <- c(1, 0)
		# warning("hiw: found a line segment that only crosses the boundary in one place.  Ignoring this segment.")

	    } 

	} # end of for 'i' loop.

	mat2 <- mat2[ -extra.points, ]
	if(all(dim(mat2) > 0)) res <- as.psp(mat2, window = res$window)

	attr(res, "side.factors") <- outvec

	return(res)

    } # end of internal 'ifun' function.

    if(verbose) {

	cat("Finding where lines from centroid to boundary along specified angles intersect with the boundary.\n")
	cat("Observed features.\n")

    }
    oseg <- apply(matrix(1:length(obs), ncol = 1), 1, ifun, x = opoly, xmid = ocen[,1], ymid = ocen[,2], l = ol, a = A, vb = verbose)
    if(verbose) cat("\nForecast features.\n")
    mseg <- apply(matrix(1:length(mod), ncol = 1), 1, ifun, x = mpoly, xmid = mcen[,1], ymid = mcen[,2], l = ml, a = A, vb = verbose)
    if(verbose) cat("Intersection points found.\n")

    thfun <- function(x) {

	a <- attributes(x)$side.factors
	id <- a != 0 & !is.na(a)

	out <- numeric(length(a)) + NA

	A <- angles.psp(x)
	out[ id ] <- sort(A)

	return(out)

    } # end of internal 'thfun' function.

    if(verbose) cat("Doing some bookeeping.\n")
    otheta <- lapply(oseg, thfun)
    mtheta <- lapply(mseg, thfun)
    if(verbose) cat("Books are in order.\n")

    # For the lengths, if the feature is not convex, and a line segment crosses in two places on the
    # same side of the centroid, need to make one of the radial segments (the shorter one) negative.

    lfun <- function(X) {

	a <- attributes(X)$side.factors
	id <- a != 0 & !is.na(a)
	A <- angles.psp(X)
	o <- order(A)
	out <- numeric(length(a)) + NA
	out[ id ] <- a[ id ] * lengths.psp(X)[ o ]

	return(out)

    } # end of internal 'lfun' function.

    if(verbose) cat("Calculating the lengths of every line segment.\n")
    obsr <- lapply(oseg, lfun)
    modr <- lapply(mseg, lfun)
    if(verbose) cat("Lengths have been calculated.\n")

    # Also need the intensities.  This is the easy part.
    intfun <- function(id, x, y) {

	ind <- as.matrix.owin(y[[ id ]])

	look <- c(x[ ind ])

	out <- c(mean(look, na.rm = TRUE), range(look, finite = TRUE))
 	names(out) <- c("mean", "min", "max")

	return(out)

    } # end of internal 'intfun' function.

    if(verbose) cat("Finding the intensity information (mean, min and max).\n")
    oint <- apply(matrix(1:length(obs), ncol = 1), 1, intfun, x = x$X, y = obs)
    mint <- apply(matrix(1:length(mod), ncol = 1), 1, intfun, x = x$Xhat, y = mod)
    if(verbose) cat("Intensity information found.\n")

    out$radial.segments <- list(X = oseg, Xhat = mseg)

    out$centers <- list(X = ocen, Xhat = mcen)
    out$intensities <- list(X = t(oint), Xhat = t(mint))
    out$angles <- list(X = otheta, Xhat = mtheta)
    out$lengths <- list(X = obsr, Xhat = modr)

    attr(out, "simplify") <- simplify

    class(out) <- "hiw"
    return(out)

} # end of 'hiw' function.

distill.hiw <- function(x, ...) {

    n1 <- dim(x$centers$X)[1]
    n2 <- dim(x$centers$Xhat)[1]

    m <- length(x$lengths$X[[ 1 ]])

    xObj <- array(NA, dim = c(m, 2, n1))
    xhatObj <- array(NA, dim = c(m, 2, n2))

    for(i in 1:n1) {

	xObj[,,i] <- matrix(x$centers$X, nrow = m, ncol = 2, byrow = TRUE) +
            cbind(x$lengths$X[[ i ]] * sin(x$angles$X[[ i ]]),
	    x$lengths$X[[ i ]] * cos(x$angles$X[[ i ]]))

    } # end of for 'i' loop.

    for(i in 1:n2) {

	xhatObj[,,i] <- matrix(x$centers$Xhat, nrow = m, ncol = 2, byrow = TRUE) + 
	    cbind(x$lengths$Xhat[[ i ]] * sin(x$angles$Xhat[[ i ]]),
            x$lengths$Xhat[[ i ]] * cos(x$angles$Xhat[[ i ]]))

    } # end of for 'i' loop.

    vxObj <- array(0, dim = c(m, 2, n1 + n2))
    vxObj[,,1:n1] <- xObj
    vxObj[,,(n1 + 1):(n1 + n2)] <- xhatObj

    attr(vxObj, "field.identifier") <- c(rep("X", n1), rep("Xhat", n2))

    return(vxObj)

} # end of 'distill.hiw' function.

summary.hiw <- function(object, ..., silent = FALSE) {

    XCenters <- object$centers$X
    XhatCenters <- object$centers$Xhat

    XInten <- object$intensities$X
    XhatInten <- object$intensities$Xhat

    n1 <- dim(XCenters)[1]
    n2 <- dim(XhatCenters)[1]

    X <- cbind(XCenters, XInten)
    Xhat <- cbind(XhatCenters, XhatInten)

    colnames(X) <- c("center.x", "center.y", "mean", "min", "max")
    colnames(Xhat) <- c("center.x", "center.y", "mean", "min", "max")

    ind <- cbind(rep(1:n1, each = n2), rep(1:n2, n1))
    colnames(ind) <- c("observed feature", "forecast feature")

    ofun <- function(id, x1, x2) {

	y1 <- x1[id[1],]
	y2 <- x2[id[2],]

	SSloc <- (y1[1] - y2[1])^2 + (y1[2] - y2[2])^2

	SSint <- (y1[3:5] - y2[3:5])^2

	res <- c(SSloc, SSint)

	names(res) <- c("SSloc", "SSavg", "SSmin", "SSmax")

	return(res)

    } # end of 'ofun' internal function.

    out <- list()

    res <- apply(ind, 1, ofun, x1 = X, x2 = Xhat)

    if(!silent) print(res)

    out$X <- X
    out$Xhat <- Xhat

    out$SS <- res

    out$ind <- ind

    invisible(out)

} # end of 'summary.hiw' function.

# summary.hiw <- function(object, ..., weights = rep(1, 7) / 7, silent = FALSE) {
# 
#     X.centers <- object$centers$X
#     Xhat.centers <- object$centers$Xhat
# 
#     n <- dim(X.centers)[1]
#     m <- dim(Xhat.centers)[1]
# 
#     if(is.null(X.centers) || is.null(Xhat.centers) || is.null(dim(X.centers)) || is.null(dim(Xhat.centers))) return(invisible())
# 
#     X.radii <- matrix(c(unlist(object$lengths$X)), nrow = n, byrow = TRUE)
#     Xhat.radii <- matrix(c(unlist(object$lengths$Xhat)), nrow = m, byrow = TRUE)
# 
#     A <- dim(X.radii)[2] # number of angles.
# 
#     X.theta <- matrix(c(unlist(object$angles$X)), nrow = n, byrow = TRUE)
#     Xhat.theta <- matrix(c(unlist(object$angles$Xhat)), nrow = m, byrow = TRUE)
# 
#     X.x <- matrix(X.centers[,1], n, A) + X.radii * sin(X.theta)
#     X.y <- matrix(X.centers[,2], n, A) + X.radii * cos(X.theta) 
#     Xhat.x <- matrix(Xhat.centers[,1], m, A) + Xhat.radii * sin(Xhat.theta)
#     Xhat.y <- matrix(Xhat.centers[,2], m, A) + Xhat.radii * cos(Xhat.theta)
# 
#     Z <- list(x = X.x, y = X.y)
#     Zhat <- list(x = Xhat.x, y = Xhat.y)
# 
#     Z <- Z$x + 1i * Z$y
#     Zhat <- Zhat$x + 1i * Zhat$y
# 
#     Zbar <- apply(Z, 1, mean, na.rm = TRUE)
#     Zhat.bar <- apply(Zhat, 1, mean, na.rm = TRUE)
# 
#     procfun <- function(ind, x1, x2, x1.int, x2.int, x1bar, x2bar, w) {
# 
# 	# z1 <- cbind(c(x1$x[ ind[1], ]), c(x1$y[ ind[1], ]))
# 	# z2 <- cbind(c(x2$x[ ind[2], ]), c(x2$y[ ind[2], ]))
# 	z1 <- x1[ ind[1], ]
# 	z2 <- x2[ ind[2], ]
# 
# 	n1 <- dim(z1)[1]
# 	n2 <- dim(z2)[1]
# 
# 	z1[ is.na(z1) ] <- 0
# 	z2[ is.na(z2) ] <- 0
# 
# 	z1.int <- x1.int[ ind[1], ]
# 	z2.int <- x2.int[ ind[2], ]
# 
# 	# z1bar <- matrix(x1bar[ ind[1], ], n1, 2, byrow = TRUE)
# 	# z2bar <- matrix(x2bar[ ind[2], ], n2, 2, byrow = TRUE)
# 	z1bar <- x1bar[ ind[1] ]
# 	z2bar <- x2bar[ ind[2] ]
# 
# 	z2star <- Conj(z2 - z2bar)
# 
# 	rhat <- abs(t(z2star) %*% (z1 - z1bar)) / ( t(z2star) %*% (z2 - z2bar) )
# 
# 	phihat <- acos(Re(t(z1 - z1bar) %*% (z2 - z2bar)) / ((t(z1 - z1bar) %*% (z1 - z1bar)) * (t(z2 - z2bar) %*% (z2 - z2bar))))
# 
# 	bhat <- z1bar - rhat * exp(-1i * phihat) * z2bar
# 
# 	zhat <- bhat + rhat * exp(-1i * phihat) * z2
# 
# 	RSS <- Re(t(Conj(z1 - zhat)) %*% (z1 - zhat))
# 
# 	SSavg <- Re(sum((z1.int[1] - z2.int[1])^2, na.rm = TRUE))
# 
# 	SSmin <- Re(sum((z1.int[2] - z2.int[2])^2, na.rm = TRUE))
# 
# 	SSmax <- Re(sum((z1.int[3] - z2.int[3])^2, na.rm = TRUE))
# 
# 	SSloc <- abs(Re(sum(bhat^2, na.rm = TRUE)))
# 
# 	SSrot <- Re(sum(phihat^2, na.rm = TRUE))
# 
# 	SSscale <- Re(sum(rhat^2, na.rm = TRUE))
# 
# 	x1trans <- z1 - z1bar
# 	x2trans <- z2 - z2bar
# 
# 	s1 <- sqrt(mean(c(x1trans)^2, na.rm = TRUE) / 2)
# 	s2 <- sqrt(mean(c(x2trans)^2, na.rm = TRUE) / 2)
# 
# 	x1scale <- x1trans / s1
# 	x2scale <- x2trans / s2
# 
# 	# theta <- atan2(sum(x1scale[,2] * x2scale[,1] - x1scale[,1] * x2scale[,2], na.rm = TRUE),
#         #                 sum(x1scale[,1] * x2scale[,1] + x1scale[,2] * x2scale[,2], na.rm = TRUE))
# 
# 	## theta <- atan2(sum(x1scale[,1] * x2scale[,1] - x1scale[,2] * x2scale[,2], na.rm = TRUE),
# 	## 		sum(x1scale[,1] * x2scale[,1] + x1scale[,2] * x2scale[,2], na.rm = TRUE))
# 
# 	# u2 <- cbind(cos(theta) * x2scale[,1] - sin(theta) * x2scale[,2], sin(theta) * x2scale[,1] + cos(theta) * x2scale[,2])
# 
# 	# RSS <- sum((u2 - x1scale)^2, na.rm = TRUE)
# 
# 	# SSscale <- mean((sqrt((x1scale[,1] - x1scale[,2])^2) - sqrt((u2[,1] - u2[,2])^2))^2, na.rm = TRUE)
# 
# 	# SSloc <- mean((colMeans(x1scale, na.rm = TRUE) - colMeans(u2, na.rm = TRUE))^2)
# 
# 	# SSrot <- theta^2
# 
# 	SStot <- w[1] * sqrt(RSS) + w[2] * SSavg + w[3] * SSmax + w[4] * SSmin + w[5] * 100 * (1 - SSscale) + w[6] * 100 * SSrot + w[7] * sqrt(SSloc)
# 
# 	res <- c(SStot, RSS, SSavg, SSmax, SSmin, SSscale, SSrot, SSloc)
# 	names(res) <- c("D", "RSS", "SSavg", "SSmax", "SSmin", "SSscale", "SSrot", "SSloc")
# 
# 	return(res)
# 
#     } # end of internal 'procfun' function.
# 
#     id <- cbind(rep(1:n, each = m), rep(1:m, n))
#     colnames(id) <- c("Obs Feature", "Fcst Feature")
# 
#     out <- t(apply(id, 1, procfun, x1 = Z, x2 = Zhat,
# 		x1.int = object$intensities$X,
# 		x2.int = object$intensities$Xhat,
# 		x1bar = X.centers, x2bar = Xhat.centers, w = weights))
# 
#     o <- order(out[,1])
# 
#     out <- cbind(id, out)
#     out <- out[o, ]
# 
#     if(!silent) print(out)
# 
#     invisible(out)
# 
# } # end of 'summary.hiw' function.

plot.hiw <- function(x, ..., which = c("X", "Xhat"), ftr.num = 1, zoom = TRUE, seg.col = "darkblue") {

    which <- match.arg(which)
    yseg <- x$radial.segments[[ which ]][[ ftr.num ]]
    if(which == "X") y <- x$X.feats[[ ftr.num ]]
    else y <- x$Y.feats[[ ftr.num ]]

    simplify <- attributes(x)$simplify

    if(is.na(simplify) || simplify == 0) ypoly <- as.polygonal(y)
    else ypoly <- as.polygonal(simplify.owin(y, dmin = simplify))

    mlab <- paste(which, " (feature ", ftr.num, ")", sep = "")

    if(zoom) {

	plot(ypoly, main = mlab, ...)
	plot(y, add = TRUE, col = "lightgray")
	plot(ypoly, add = TRUE)

    } else {

	plot(y, main = mlab, ...)
	plot(ypoly, add = TRUE)

    } # end of if else 'zoom' stmts.

    plot(yseg, add = TRUE, col = seg.col)

    invisible()

} # end of 'plot.hiw' function.

print.hiw <- function(x, ...) {

    print(x$identifier.function)
    print(x$identifier.label)

    print(x$radial.segments)

    print(x$centers)

    print(x$intensities)

    print(x$angles)

    print(x$lengths)

    invisible()

} # end of 'print.hiw' function.

# helmerter <- function(n) {
# 
#     hfun <- function(x, n) {
# 
# 	hr <- -1 / sqrt(x * (x + 1))
# 	if(x < n) res <- c(rep(hr, x - 1), -x * hr, rep(0, n - x))
# 	else res <- c(rep(hr, x - 1), -x * hr)
# 	return(res)
# 
#     } # end of internal 'hfun' function.
# 
#     r <- matrix(2:n, ncol = 1)
#     out <- t(apply(r, 1, hfun, n = n))
#     out <- rbind(rep(1 / sqrt(n), n), out)
#     return(out)
# 
# } # end of 'helmerter' function.

# preshaper <- function(centers, radii, angles) {
# 
#     x <- centers[,1] + radii * sin(angles)
#     y <- centers[,2] + radii * cos(angles)
# 
#     z <- x + 1i * y
# 
#     k <- length(z)
#     H <- helmerter(n = k)[-1,]
# 
# 
# } # end of 'preshaper' function.
