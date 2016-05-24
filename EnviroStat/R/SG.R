# These R functions are based on the Splus codes provided by Prof. Paul Sampson.

# No changes from Version 0.0


Falternate3 <- function(disp, coords, model = 1., a0 = 0.1, t0 = 0.5, 
     max.iter = 50., max.fcal= 100., alter.lim = 50., tol = 1e-05, prt = 0., 
     dims = 2., lambda = 0., ncoords, dev.mon = NULL, verbose = FALSE)
{
	# Do simultaneous estimation of coords and exponential or gaussian 
	# variogram by alternating weighted least squares.
	# This version permits dimension > 2 for scaling.
	# In the plotting we'll use a plot symbol proportional to the
	# third coordinate.
	#
	# This version also passes a smoothing parameter to the optimization.
	# This parameter probably is not scaled exactly the same as it is
	# in sinterp and this has not been investigated yet.  
	# Warning: make sure that coords are scaled reasonably small
	# before attempting to compute; otherwise matrix inversion is
	# likely not to work in calculation of bending energy matrix.
	#
	# Other arguments:
	#       model:   1 (exponential), 2 (gaussian)
	#       a0,t0:   initial variogram parameter estimates
	#       max.iter, max.fcal:  control parameters for calls to nlmin
	#                (same values used in MDS step and in variogram step)
	#       dev.mon: device number for plot monitoring convergence of objective
	#       ncoords: optional initial coordinates to use (2 dim) if not G-plane
	#
    if (!is.null(dev.mon)) {
        dev.plot <- dev.cur()
        if (isTRUE(dev.mon)) {
            dev.mon <- getOption("device")
        }
        dev.fn <- match.fun(dev.mon)
        dev.fn()
        dev.mon <- dev.cur()
        dev.set(dev.plot)
    }
    maxdiff = 1
	par(mfrow = c(1., 2.))
	if(missing(ncoords))
		ncoords <- coords
	# Initialize 3rd coord with small random numbers.
	# Fix z coords of first three points to be zero.
	BM <- Fbenergy(coords)
	if(dims > 2.) {
		ncoords <- cbind(coords, matrix(rnorm(coords * (dims - 2.),
			0., sqrt(var(coords[, 1.])/100.)), nrow = nrow(coords),
			ncol = (dims - 2.)))
		ncoords[1.:3., 3.] <- 0.
	}
    plotCoordChange(coords, ncoords, dims, 0, maxdiff, lambda)

	BEpenalty <- lambda * sum(diag(t(ncoords) %*% BM %*% ncoords))
	h <- Fdist(ncoords)
	h.lt <- h[row(h) < col(h)]
	#
	# Initial fit
	disp.lt <- disp[row(disp) < col(disp)]
	i <- 0.
	variogfit <- Fvariogfit3(disp.lt, h.lt, model = model, a0 = a0, t0 = t0,
		bep = BEpenalty, max.iter = max.iter, max.fcal = max.fcal, verbose = verbose)
    objf <- variogfit$objf

    if (!is.null(dev.mon)) {
        ## Set up plot for monitoring optimization criterion
        dev.set(dev.mon)
        plot(c(0., alter.lim), c(0., objf), type = "n", xlab = "iteration",
             ylab = "penalized root mean square relative error")
        points(0., objf)
        dev.set(dev.plot)
    }
	a0 <- variogfit$a[1.]
	t0 <- variogfit$t0
    plotDispersionDist(h.lt, disp.lt, model, variogfit, ncoords)

	maxdiff = 100.
	while(i <= alter.lim & maxdiff > tol) {
		#
		i <- i + 1.
		oldcoords <- ncoords
		ncoords <- Fmdsfit3(disp.lt, ncoords, model = model, a = 
			variogfit$a, t0 = variogfit$t0, prt = prt, dims = dims,
			lambda = lambda, bem = BM, max.iter = max.iter, 
			max.fcal = max.fcal, verbose = verbose)
		objf <- ncoords$objf
        if (!is.null(dev.mon)) {
            dev.set(dev.mon)
            text(i, objf, "m")
            dev.set(dev.plot)
        }
		ncoords <- ncoords$ncoords
		maxdiff <- max(abs(ncoords - oldcoords))

        plotCoordChange(coords, ncoords, dims, i, maxdiff, lambda)

		BEpenalty <- lambda * sum(diag(t(ncoords) %*% BM %*% ncoords))
		h <- Fdist(ncoords)
		h.lt <- h[row(h) < col(h)]
		variogfit <- Fvariogfit3(disp.lt, h.lt, model = model, a0 = a0,
			t0 = t0, bep = BEpenalty, max.iter = max.iter, max.fcal
			 = max.fcal)
		objf <- variogfit$objf
        if (!is.null(dev.mon)) {
            dev.set(dev.mon)
            text(i, objf, "v")
            dev.set(dev.plot)
        }
        a0 <- variogfit$a[1.]
        t0 <- variogfit$t0
        plotDispersionDist(h.lt, disp.lt, model, variogfit, ncoords)
	}
	list(variogfit=variogfit, ncoords=ncoords)
}


plotCoordChange <- function(coords, ncoords, dims, i, maxdiff, lambda) {
	temp <- setplot(rbind(coords, ncoords[, 1.:2.]), axes = T)
	text(coords,labels=c(1:(dim(coords)[1])))
	if(dims == 3.) {
		symbols(ncoords[, 1.], ncoords[, 2.],
                circles = (ncoords[,3.] - min(ncoords[, 3.]) + 0.01),
                inches = 0.15, add = T)
	}
    suppressWarnings(
        arrows(coords[, 1.], coords[, 2.],
               ncoords[, 1.], ncoords[, 2.], length=0.075) )
    title(main = paste("D-plane, iteration ", i),
          sub = paste("Max coordinate change =", format(round(maxdiff, 4.))) ,
          xlab=substitute(list(lambda) == list(x), list(x=lambda)) )
	par(pin = temp$oldpin)
}


plotDispersionDist <- function(h.lt, disp.lt, model, variogfit, ncoords) {
    rmsre <- sqrt(mean(((disp.lt - variogfit$fit)/variogfit$fit)^2))
    plot(h.lt, disp.lt, xlab = "D-plane distance",
         ylab = paste("Dispersion, rmsre =", format(round(rmsre, 4.))),
         xlim = c(0., max(h.lt)), ylim = c(0., max(disp.lt)))
    title(main = paste(ifelse(model == 1., " Exponential",
                              " Gaussian"), "Fitted Variogram"), cex = 0.75)
    title(main = paste("\ncriterion = ", format(round(variogfit$objf, 4.))),
          cex = 0.75)
    abline(h = 2., lty = 2.)
    
    lines(h.lt[order(h.lt)], variogfit$fit[order(h.lt)])
    
    points(0., variogfit$a[1.], pch = 1.)
}

Ftransdraw <- function(disp, Gcrds, MDScrds, gridstr, sta.names, lambda = 0., lsq = FALSE, eye,
	model = 1., a0 = 0.1, t0 = 0.5)
{
# Purpose: for varying (user-supplied) values of the spline smoothing
#          parameter lambda, 
# - compute image of coordinates, interpoint distances, and dist-disp plot.
# - compute and draw image of regular grid (2D or 3D perspective plot)
       # Computed prior to execution of this code:
       # - disp: spatial dispersion matrix
       # - Gcrds: geographic coordinates (nx2)
       # - MDScrds:  kyst mds solution (nx2 or nx3) - using Falternate3
       # - gridstr: regular grid structure on the G-plane (from Fmgrid)
	
	
	on.exit(par(par0))
	par0 <- par(no.readonly = TRUE)
	oldpar <- par(xpd = T)
	grid <- gridstr$grid
	if(missing(sta.names)) 
            sta.names <- format(1.:nrow(Gcrds))
	rl.ind <- gridstr$rl.ind
	Ddim <- ncol(MDScrds)
	#
	grid.nm <- !#
	is.na(grid[, 1.])
	plim <- setplot(grid[, 1.], grid[, 2.], axes = T)
	text(Gcrds[, 1.], Gcrds[, 2.], sta.names, cex = 0.75)
	lines(grid[, 1.], grid[, 2.])
	title(main = "Geographic Coordinates")
	cat("Click anywhere on plot to continue (Left button)\n")
	junk <- locator(n=1)
	#
	par(pin = plim$oldpin)
	if(Ddim == 2.)
		par(mfrow = c(1., 2.))
	while(length(lambda) == 1.) {
		Tspline <- sinterp(t(Gcrds), MDScrds, lam = lambda, lsq = lsq)
		Dcrds <- t(seval(t(Gcrds), Tspline)$y)
		Ddist <- Fdist(Dcrds)
		Dgrid <- matrix(NA, nrow(grid), Ddim)
		Dgrid[grid.nm,  ] <- t(seval(t(grid[grid.nm,  ]), Tspline)$y)
		h.lt <- Ddist[row(Ddist) < col(Ddist)]
		disp.lt <- disp[row(disp) < col(disp)]
		variogfit <- Fvariogfit3(disp.lt, h.lt, model = model, a0 = a0,
			t0 = t0)
		a0 <- variogfit$a[1.]
		t0 <- variogfit$t0
		rmsre <- round(sqrt(mean(((disp.lt[order(h.lt)] - variogfit$
			fit[order(h.lt)])/variogfit$fit[order(h.lt)])^2.)),
			4.)
		plot(h.lt, disp.lt, xlab = paste("D-space distance"), 
                   ylab = paste("Dispersion, rmsre =", format(
			rmsre)), ylim = c(0., max(disp.lt)), xlim = c(0., max(
			h.lt)))
		title(main = paste(ifelse(model == 1.,
			" Exponential", " Gaussian"), "Variogram"), cex = 0.75)
                if(max(disp.lt)>2) abline(h=2., lty=2.)
		lines(h.lt[order(h.lt)], variogfit$fit[order(h.lt)])
		points(0., variogfit$a[1.], pch = 1.)
                title(sub = substitute(list(lambda) == list(x), list(x=lambda))  )
	#	cat("Click anywhere on plot to continue (left button)   \n")
	#	junk <- locator(n=1)
		if(Ddim == 2.) {
		plim <- setplot(Dgrid[, 1.], Dgrid[, 2.], axes = T)
		text(Dcrds[, 1.], Dcrds[, 2.], sta.names, cex = 0.75)
		lines(Dgrid[, 1.], Dgrid[, 2.])
		title(main = "D-plane Coordinates", 
                    xlab = substitute(list(lambda) == list(x), list(x=lambda)),
                sub = expression(paste("[",lambda," : Spline smoothing parameter]") ))
		       par(pin = plim$oldpin)
		}
		else {
            stop("3-D perspective plot is not yet implemented")
#			Dxrange <- range(Dgrid[grid.nm, 1.])
#			Dyrange <- range(Dgrid[grid.nm, 2.])
#			if(missing(eye))
#			eye <- c(Dxrange[1.] - 4. * (Dxrange[2.] - 
#				Dxrange[1.]), Dyrange[1.] - 4. * (Dyrange[2.] - Dyrange[1.]), 
#                                    max(Dgrid[grid.nm, 3.]) + 5. * (max(Dgrid[grid.nm,
#					3.]) - min(Dgrid[grid.nm, 3.])))
#			while(length(eye) == 3.) {
#				Dinterp <- interp(Dgrid[grid.nm & rl.ind == 1.,	1.], 
#                                Dgrid[grid.nm & rl.ind == 1.,2.], 
#                                Dgrid[grid.nm & rl.ind == 1.,3.], 
#                                xo = seq(Dxrange[1.], Dxrange[2.], length = 10.), 
#                                yo = seq(Dyrange[1.], Dyrange[2.], length = 10.))
#				Dinterp$z[is.na(Dinterp$z)] <- 0.
#				lp.out <- persp(Dinterp$x, Dinterp$y, Dinterp$
#					z, eye = eye,
#                    lty=1:3, col=rep(1L, 3), lwd=rep(0L, 3))
#				Dcrds.pp <- trans3d(Dcrds[, 1.], Dcrds[, 2.],
#					Dcrds[, 3.], lp.out)
#				#
#				text(Dcrds.pp$x, Dcrds.pp$y, sta.names, cex = 
#					0.75)
#				Dgrid.pp <- Dgrid[, 1.:2.]
#				temp.pp <- trans3d(Dgrid[grid.nm, 1.], 
#                                              Dgrid[grid.nm, 2.], 
#                                              Dgrid[grid.nm, 3.], lp.out)
#				Dgrid.pp[grid.nm,  ] <- cbind(temp.pp$x, temp.pp$y)
#				lines(Dgrid.pp[, 1.], Dgrid.pp[, 2.])
#				title(main = paste(
#					"D-space Coordinates,  eye = ", 
#                                           round(eye[1.]), round(eye[2.]), 
#                                           round(eye[3.]), sep = " "), 
#                    xlab = substitute(list(lambda) == list(x), list(x=lambda)),
#                    sub = expression(paste("[",lambda," : Spline smoothing parameter]") ))
#		  
#				cat("Click anywhere on plot to continue\n")
#				junk <- locator(n=1)
#				cat("Enter 3 numbers for new eye perspective\n"
#					)
#				eye <- scan(n = 3.)  
#			}
		}
             cat("Enter value for new lambda (Hit return to stop) \n")
	     lambda <- scan(n=1)  
             if (length(lambda) >0) {
                cat("Click anywhere on plot to continue (left button)  \n")
		junk <- locator(n=1) 
                    }  
}
	list(Dcrds=Dcrds, Ddist=Ddist)
}

Fbenergy <- function(crds)
{
	# Function to compute the bending energy matrix for thin-plate spline
	# mappings of the coordinates in the nx2 matrix crds.
	# Follow notation of Mardia, Kent & Walder, 1991.
	n <- nrow(crds)
	Dist2 <- Fdist(crds)^2.
	Sig <- Dist2 * logb(Dist2)
	Sig[row(Sig) == col(Sig)] <- 0.
	Tc <- t(cbind(1., crds))
	P <- t(Tc) %*% solve(Tc %*% t(Tc)) %*% Tc
	IP <- diag(n) - P
	B1 <- IP %*% Sig %*% IP
	B <- FMPinv(B1)
	return(B)
}


Fdist <- function(crds)
{
	# Function to compute interpoint distances for nxp coordinate matrix.
	dist <- matrix(0., nrow(crds), nrow(crds))
	for(i in 1.:ncol(crds)) {
		dist <- dist + outer(crds[, i], crds[, i], "-")^2.
	}
	dist <- sqrt(dist)
	return(dist)
}

Feiggrid <- function(grid, xmat, coef)
{
	#
	#	Function to compute and return eigenvalues/vectors of affine
	#	derivative matrix at locations in grid for spline given by
	#	coefficients in coef.
	#
	if(!is.matrix(grid)) stop(paste(substitute(grid), "is not a matrix"))
	if(!is.matrix(xmat))
		stop(paste(substitute(xmat), "is not a matrix"))
	if(!is.matrix(coef))
		stop(paste(substitute(coef), "is not a matrix"))
	datain <- matrix(as.single(xmat), nrow(xmat), ncol(xmat))
	ndat3 <- nrow(coef)
	nr <- ncol(coef)
	coefin <- matrix(as.single(coef), ndat3, nr)
	igrid <- !is.na(grid[, 1.])
	ngrid <- nrow(grid[igrid,  ])
	cgrid <- matrix(as.single(grid[igrid,  ]), ngrid, ncol(grid))
	eigvec <- matrix(as.single(cbind(cgrid, cgrid)), nrow = ngrid)
	eigval <- matrix(as.single(cgrid), nrow = ngrid)
	Q <- .Fortran("eiggrid",
		ngrid,
		cgrid,
		ndat3,
		nr,
		datain,
		coefin,
		eigvec = eigvec,
		eigval = eigval)
	eigvec <- cbind(grid, grid)
	eigvec[igrid,  ] <- Q$eigvec
	eigval <- grid
	eigval[igrid,  ] <- Q$eigval
	list(grid=grid, eigvec=eigvec, eigval=eigval)
}


Flamb2 <- function(geoconfig, latrf1 = NA, latrf2 = NA, latref = NA, lngref = NA)
{
 # Evaluate Lambert projection for geoconfig: (lat, -long)

	geo <- as.matrix(geoconfig)
	if(dim(geo)[[2.]] != 2.)
		stop("the input should be an nx2 matrix")
	xy.coord <- array(0., dim(geo))
	if(is.na(latref))
		latref <- sum(range(geo[, 1.]))/2.
	if(is.na(lngref))
		lngref <- sum(range(geo[, 2.]))/2.
	if(is.na(latrf1))
		latrf1 <- latref - (3./10.) * diff(range(geo[, 1.]))
	if(is.na(latrf2))
		latrf2 <- latref + (3./10.) * diff(range(geo[, 1.]))
	lat <- geo[, 1.]
	long <- geo[, 2.]
	pi <- 3.14159265
	a <- 6378137.
	b <- 6356752.
	radlf1 <- (pi/180.) * latrf1
	radlf2 <- (pi/180.) * latrf2
	radlgf <-  - (pi/180.) * lngref
	radltf <- (pi/180.) * latref
	eccen <- sqrt((a^2. - b^2.)/a^2.)
	capr <- (a * (1. - eccen^2.))/((1. - eccen^2. * sin(radltf)^2.)^1.5)
	n <- logb(cos(radlf1)/cos(radlf2))/(logb(tan(pi/4. + radlf2/2.)/tan(
		pi/4. + radlf1/2.)))
	capf <- (cos(radlf1) * ((tan(pi/4. + radlf1/2.))^n))/n

	rho0 <- (capr * capf)/((tan(pi/4. + radltf/2.))^n)
	radlat <- (pi/180.) * lat
	radlng <-  - (pi/180.) * long
	theta <- n * (radlng - radlgf)
	rho <- (capr * capf)/((tan(pi/4. + radlat/2.))^n)
	x <- (0.001) * rho * sin(theta)
	y <- (0.001) * (rho0 - rho * cos(theta))
	xy <- cbind(x, y)
	list(xy =xy, latrf1=latrf1, latrf2=latrf2, latref=latref, lngref=lngref)

}


Fmdsfit3 <- function(disp.lt, coords, model = 1., a, t0, max.iter = 25., max.fcal = 100.,
	prt = 0., dims = 2., lambda = 0., bem, verbose = FALSE)
{
	# Fit a new set of coordinates for given variogram parameters
	# Data to be made available to the least squares program
	# are now the variogram parameters (model=1: exp, model=2: gauss)
	# and dispersions.

        L=lambda
        BM= bem
        Y = disp.lt
        MODEL = list(T0 = t0, A = a, M = model)
        DIM= dims
        X0 = t(coords[1.:2.,  ])
	n <- nrow(coords)
	if(dims == 3.) {
		x <- c(coords[3., -3.], t(coords[4.:n,  ]))
	}
	else {
		x <- c(t(coords[ - c(1., 2.),  ]))
	}
    tt <- nlm(Tressq.mds3, x, L=L, BM=BM,MODEL = MODEL, Y=Y, DIM=DIM, X0=X0)
    temp <- list(x=0, converged = FALSE, conv.type = 1)
    temp$x <- tt$estimate
    if (tt$code<=2) temp$converged <- T
    temp$conv.type <- tt$code
	#
	if (verbose) cat("MDS-convergence: ", temp$converged, "\n")
	#
	if (verbose) cat("  nlm code =  ", temp$conv.type, "\n")
	objf <- Tressq.mds3(temp$x, L=L, BM=BM,MODEL = MODEL, Y=Y, DIM=DIM, X0=X0)
	# Reassemble
	#
	if (verbose) cat("   criterion: ", objf, "\n")
	if(dims == 2.) {
		ncoords <- matrix(c(t(coords[1.:2.,  ]), temp$x), byrow = T,
			ncol = 2.)
	}
	if(dims == 3.) {
		ncoords <- matrix(c(t(coords[1.:2.,  ]), temp$x[1.:2.], 0.,
			temp$x[ - (1.:2.)]), byrow = T, ncol = 3.)
	}
	list(objf=objf, ncoords=ncoords)
}


Fmgrid <- function(xlim, ylim, xn = 8., xspace, xres, yn = 8., yspace, yres)
{
	# Function to generate points on a grid.  Points are assembled
	# in an nx2 matrix with NA's separating series of verticle
	# and horizontal lines in the grid.
	# - xn and yn specify the number of vertical and horizontal lines, 
	# respectively.  These parameters are overridden by xspace and yspace
	# if specified.
	# - xspace and yspace specify the distance between successive
	# vertical and horizontal lines, respectively.
	# - xres and yres specify the distance between points generated
	# along horizontal and vertical lines, respectively.
	# - if xres and yres are not specified, then points are generated
	# only at the nodes of intersection of the vertical and horizontal
	# lines.  Note that these nodes appear in duplicate as sequences
	# of points are generated first for the vertical lines and then
	# for the horizontal lines.
	if(length(xlim) < 2.) stop(paste(substitute(xlim), 
			"is not a 2 element vector"))
	if(length(ylim) < 2.)
		stop(paste(substitute(ylim), "is not a 2 element vector"))
	xrange <- xlim[2.] - xlim[1.]
	yrange <- ylim[2.] - ylim[1.]
	if(missing(xspace))
		xspace <- xrange/xn
	if(missing(yspace))
		yspace <- yrange/yn
	if(missing(xres))
		xres <- xspace
	if(missing(yres))
		yres <- yspace
	xl <- c(seq(xlim[1.], xlim[2.], xspace), NA)
	xd <- c(seq(xlim[1.], xlim[2.], xres), NA)
	yl <- c(seq(ylim[1.], ylim[2.], yspace), NA)
	yd <- c(seq(ylim[1.], ylim[2.], yres), NA)
	nxl <- length(xl)
	nxd <- length(xd)
	nyl <- length(yl)
	nyd <- length(yd)

   xc <- c(rep(xl[1:(nxl-1)], each=nyl), rep(xd, nyl - 1.))
   yc <- c(rep(yd, nxl - 1.), rep(yl[1:(nyl-1)], each=nxd))
	rl.ind <- c(rep(1., (nxl - 1.) * nyd), rep(2., nxd * (nyl - 1.)))
	xc[is.na(yc)] <- NA
	yc[is.na(xc)] <- NA
	grid <- cbind(xc, yc)
	list(grid=grid, rl.ind=rl.ind)
}
FMPinv <- function(X)
{
	# Generalized inverse of symmetric X
	Xeig <- eigen(X, symmetric = T)
	Xval <- Xeig$values
	ipos <- Xval > 1e-08
	if(sum(ipos) > 1.) {
		X1 <- Xeig$vectors[, ipos] %*% diag(1./Xval[ipos]) %*% t(Xeig$
			vectors[, ipos])
	}
	else {
		if(sum(ipos) == 1.) {
			X1 <- Xeig$vectors[, ipos] %*% matrix(Xeig$vectors[
				, ipos], nrow = 1.)/Xval[ipos]
		}
		else {
			X1 <- matrix(0., nrow(X), ncol(X))
		}
	}
	Xeig <- eigen(X1, symmetric = T)
	Xsel <- abs(Xeig$values) < 100000000.
	X1 <- Xeig$vectors[, Xsel] %*% diag(Xeig$values[Xsel]) %*% t(Xeig$
		vectors[, Xsel])
	return(X1)
}



Fvariogfit3 <- function(disp.lt, h.lt, model = 1., a0 = 0.1, t0 = .5,
                 max.iter = 25., max.fcal= 100., bep = 0., verbose = FALSE)
{
	# Fit an exponential or gaussian variogram
 
        Y = disp.lt
        H = h.lt
        if (t0 > 700) t0 = .1
        if (a0 > 2) a0 = .1
        MODEL = list(M=model)
  
        BEP = bep        

    tt <- nlm(Tressq.variog3, logb(c(t0, a0)), MODEL=MODEL, H=H, Y=Y, BEP=bep)
    temp <- list(x=0, converged = FALSE, conv.type = 1)
    temp$x <- tt$estimate
    if (tt$code<=2) temp$converged <- T
    temp$conv.type <- tt$code
	t0 <- exp(temp$x[1.])
	a0 <- exp(temp$x[2.])
	a1 <- 2. - a0
	if (verbose) cat("VAR-convergence: ", temp$converged, "\n")
	if (verbose) cat("  nlm code =  ", temp$conv.type, "\n")
	objf <- Tressq.variog3(temp$x, MODEL=MODEL, H=H, Y= Y, BEP=bep)
	if (verbose) cat("    criterion: ", objf, "\n")
	if(model == 1.) {
		fit <- a0 + a1 * (1. - exp( - t0 * h.lt))
	}
	else {
		fit <- a0 + a1 * (1. - exp( - t0 * h.lt^2.))
	}
	a <- c(a0,a1)
	list(objf=objf, t0=t0, a=a , fit=fit)
}

integ <- function(start, xmat, coef, ind = 2., xlimt, iter.limit = 10.)
{
	if(length(start) != 2.)
		stop(paste(substitute(start), "is not a 2 element vector"))
	if(!is.matrix(xmat))
		stop(paste(substitute(xmat), "is not a matrix"))
	if(!is.matrix(coef))
		stop(paste(substitute(coef), "is not a matrix"))
	if(missing(xlimt))
		xlimt <- c(apply(xmat[, 1.:2.], 2., range))
	if(length(xlimt) != 4.)
		stop(paste(substitute(xlimt), "is not a 4 element vector"))
	nn <- as.integer(1000.)
	start <- as.single(start)
	xylim <- as.single(xlimt)
	datlen <- nrow(xmat)
	datain <- matrix(as.single(xmat), datlen, ncol(xmat))
	coefdim <- ncol(coef)
	coefin <- matrix(as.single(coef), datlen + 3., coefdim)
	ind <- as.integer(ind)
	grid <- matrix(as.single(0.), 2., nn)
	fldmag <- vector("single", nn)
	ngrid <- vector("integer", 1.)
	iter.limit <- as.integer(iter.limit)
	w1 <- grid
	w2 <- fldmag
	w3 <- fldmag
	Q <- .Fortran("integ",
		nn = nn,
		start,
		xylim,
		datain,
		datlen,
		coefin,
		coefdim,
		grid = grid,
		ngrid = ngrid,
		ind,
		fldmag = fldmag,
		iter.limit,
		w1,
		w2,
		w3)
	R <- list(grid = Q$grid, ngrid = Q$ngrid, fldmag = Q$fldmag, nn = Q$nn)
	if(R$nn == 999.)
		stop("internal grid buffer overflow")
	R$grid <- R$grid[, c(1.:R$ngrid)]
	R$fldmag <- R$fldmag[c(1.:R$ngrid)]
	R$fldmag[R$grid[1.,  ] == 999.] <- NA
	R$grid[R$grid == 999.] <- NA
	return(R)
}


setplot <- function(xdata, ydata, pretty.call = TRUE, maxdim, axes = FALSE)
{
	if(missing(xdata))
		stop("no xdata nor ydata was passed to setplot")
	if(is.matrix(xdata)) {
		if(ncol(xdata) != 2.)
			stop(paste(substitute(xdata), "has too many columns"))
		ydata <- xdata[, 2.]
		xdata <- xdata[, 1.]
	}
	else if(is.list(xdata)) {
		ydata <- xdata$y
		if(is.null(ydata))
			stop(paste(substitute(xdata), "has no y component"))
		xdata <- xdata$x
		if(is.null(xdata))
			stop(paste(substitute(xdata), "has no x component"))
	}
	else if(missing(ydata))
		stop("no ydata was passed to setplot")
	if(pretty.call) {
		xdata <- pretty(xdata)
		ydata <- pretty(ydata)
	}
	xlim <- range(xdata)
	ylim <- range(ydata)
	xrng <- xlim[2.] - xlim[1.]
	yrng <- ylim[2.] - ylim[1.]
	prng <- max(xrng, yrng)
	oldpin <- par("pin")
	if(missing(maxdim)) {
		if(xrng/yrng > oldpin[1.]/oldpin[2.]) {
			maxdim <- oldpin[1.]
			prng <- xrng
		}
		else {
			maxdim <- oldpin[2.]
			prng <- yrng
		}
	}
	newpin <- (maxdim * c(xrng, yrng))/prng
	par(pty = "m", pin = newpin)
	plot(xlim, ylim, type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = 
		"i", axes = axes)
    if (interactive()) {
        warning("Plot limits changed.  Reset before making new plots!")
        warning(paste("Old limits: ", paste(oldpin, collapse = "  ")))
        warning(paste("New limits: ", paste(newpin, collapse = "  ")))
            list(xlim=xlim, ylim=ylim, oldpin=oldpin, newpin=newpin)
    }
}


seval <- function(x, tpsp)
{
	tpsp.names <- c("x", "y", "m", "lam", "lsq", "b", "sol", "ainf", "linf",
		"f", "a")
	if(!all(names(tpsp) == tpsp.names))
		stop(paste(substitute(tpsp), "is not in the correct format"))
	dblx <- tpsp$x
	d <- nrow(dblx)
	if(is.matrix(x)) {
		dp <- nrow(x)
		if(dp != d) {
			dp <- ncol(x)
			if(dp == d)
				x <- t(x)
		}
		np <- ncol(x)
	}
	else {
		dp <- length(x)
		np <- as.integer(1.)
	}
	if(dp != d)
		stop(paste(substitute(x), 
			"dimension incompatible with node matrix dimension"))
	n <- ncol(dblx)
	x <- matrix(as.double(x), d, np)
	nna.index <- apply(!is.na(x), 2., all)
	nna <- sum(nna.index)
	x.index <- matrix(nna.index, d, np, byrow = TRUE)
	x.nna <- matrix(x[x.index], d, nna)
	if((d %/% 2.) * 2. != d)
		odd <- as.integer(1.)
	else odd <- as.integer(0.)
	m <- tpsp$m
	lam <- tpsp$lam
	nlam <- as.integer(length(lam))
	sdim <- as.integer(c(dim(tpsp$sol), 1.)[1.:3.])
	sol <- array(tpsp$sol, dim = sdim)
	slen <- sdim[1.]
	ny <- sdim[2.]
	f <- tpsp$f
	flen <- as.integer(length(f))
#	dbla <- tpsp$a
 	dbla <- tpsp$a[,1]
	alen <- as.integer(length(dbla))
	ximxj <- vector("double", d)
	yout <- as.double(0.)
	ydim <- vector("integer", 3.)
	ydim[1.] <- ny
	ydim[2.] <- nna
	ydim[3.] <- nlam
	y <- array(as.double(0.), dim = ydim)
	R <- .Fortran("seval",
		x.nna,
		d=as.integer(d),
		nna=as.integer(nna),
		dblx,
		n=as.integer(n),
		odd,
		m=as.integer(m),
		lam,
		nlam=as.integer(nlam),
		sol,
		slen=as.integer(slen),
		f,
		flen,
		dbla,
		alen=as.integer(alen),
		ximxj,
		yout = yout,
		y = y,
		ny=as.integer(ny))
	if(nlam == 1.)
		ydim <- ydim[-3.]
	ydim[2.] <- np
	x.index <- array(x.index, dim = ydim)
	y <- array(as.double(NA), dim = ydim)
	y[x.index] <- R$y
	list(x=x, y=y)
}


sinterp <- function(x, y, m = 2., lam = 0., lsq = FALSE)
{
	if(!is.matrix(x))
		stop(paste(substitute(x), "is not a matrix"))
	if(!is.matrix(y))
		stop(paste(substitute(y), "is not a matrix"))
	m <- as.integer(m)
	lam <- as.double(lam)
	lsq <- as.logical(lsq)
	n <- ncol(x)
	yr <- nrow(y)
	if(n != yr) {
		n <- nrow(x)
		if(n != yr)
			stop(paste(substitute(x), "and", substitute(y), 
				"are not conformable matrices"))
		else x <- t(x)
	}
	d <- nrow(x)
	ny <- ncol(y)
	x <- matrix(as.double(x), d, n)
	y <- matrix(as.double(y), n, ny)
	nlam <- as.integer(length(lam))
	f <- vector("integer", (d + m))
	f[1.] <- as.integer(1.)
	for(i1 in 2.:(d + m))
		f[i1] <- as.integer(f[i1 - 1.] * (i1 - 1.))
	alen <- as.integer(f[d + m]/f[m]/f[d + 1.] + n)
	alenp <- as.integer((alen * (alen + 1.))/2.)
	dlen <- as.integer((alenp * nlam) + d + (n + (2. + n) * (alen - n)))
	soldim <- vector("integer", 3.)
	soldim[1.] <- alen
	soldim[2.] <- ny
	soldim[3.] <- nlam
	a <- matrix(as.double(0.), alenp, nlam)
	b <- matrix(as.double(0.), (alen - n), ny)
	sol <- array(as.double(0.), dim = soldim)
	iwrk <- vector("integer", alen)
	ainf <- vector("integer", nlam)
	linf <- vector("integer", ny)
	db <- vector("double", dlen)
	R <- .Fortran("intdrv",
		x = x,
		d = as.integer(d),
		n = as.integer(n),
		y = y,
		ny = as.integer(ny),
		m =  as.integer(m),
		alen = as.integer(alen),
		f = as.integer(f),
		a = a,
		b = b,
		sol = sol,
		iwrk = as.integer(iwrk),
		lam = lam,
		nlam = as.integer(nlam),
		ainf = ainf,
		linf = as.integer(linf),
		lsq = as.logical(lsq),
		db,
		dlen = as.integer(dlen))
	if(nlam == 1.)
		R$sol <- array(R$sol, dim = soldim[-3.])
	junk <- list(x = R$x, y = R$y, m = R$m, lam = R$lam, lsq = R$lsq, b = R$b, 
                  sol = R$sol, ainf = R$ainf, linf = R$linf, f = R$f, a = R$a)
invisible(junk)
}

Tressq.mds3 <- function(x, L, BM, Y, MODEL, DIM, X0)
{  

	# Function to evaluate weighted residual sum of squares for
	# exponential or gaussian variogram when fitting site coordinates for
	# fixed variogram parameters. 
	
	# Response (dispersion) is assumed in the vector 'Y' (lt part), 
	# Variogram parameters are in 'MODEL', including indicator 'M'
	# of exponential or gaussian variogram.
	# Bending energy matrix is in 'BM'; penalty weight in 'L'.
	
	# Reassemble
	if(DIM == 3.) x <- c(x[1.:2.], 0., x[ - c(1.:2.)])
	crds <- #
	matrix(c(X0, x), byrow = T, ncol = DIM)
	BEP <- L * sum(diag(t(crds) %*% BM %*% crds))
	dist <- Fdist(crds)
	dist.lt <- dist[row(dist) < col(dist)]
	G <- cbind(1., 1. - exp( - MODEL$T0 * dist.lt^MODEL$M))
	#
	fit <- G %*% MODEL$A
        resid <- (Y - fit)
	#
	objf <- sum((resid)^2.) + BEP
	return(objf)
}


Tressq.variog3 <- function(x, MODEL, H, Y, BEP)
 {
	# Function to evaluate weighted residual sum of squares for
	# exponential or gaussian variogram. 
	# Response (dispersion) is assumed in the vector 'Y' (lt part), 
	# Distance is assumed in vector 'H' (also lt part).
	# Model indicator is assumed in list 'MODEL'.
	# Bending energy penalty is in 'BEP'.
	
	# Parameters are constrained to be positive by representing them
	# as exp().
	t0 <- exp(x[1.])
	a0 <- # postive nugget
	exp(x[2.])
	a1 <- 2. - # asymptote (sill) of 2
	(a0)
	a <- c(a0, a1)
	if(MODEL$M == 1.) {
		G <- cbind(1., 1. - exp( - t0 * H))
	}
	else {
		G <- cbind(1., 1. - exp( - t0 * H^2.))
	}
	#
	fit <- G %*% a
        resid <- (Y - fit)
	#
	objf <- sum((resid)^2.) + BEP
	return(objf)
}

Disp.link <- function(disp.mx, coords, dmap, ddisp, names, device = getOption('device'))
 {
	# Function to plot and link points in a dispersion-distance plot and a
	# geographic map.
	# User may identify points on the dispersion scatter in order to
	# identify line segments on the coordinate plot, 
	# **and/or** identify individual sites on the coordinate plot in order
	# to highlight the corresponding set of points on the dispersion scatter.
	#
	# Uses 'setplot' to set up coordinates for geographic map
	# and 'Fdist' to compute distances from coords.
	# disp.mx should be nxn matrix of dispersions (i.e. Var(Z(x)-Z(y))
	# coords should be nx2 matrix of coordinates
	# Optional input: 
	#       dmap=device number for existing window to be used for map
	#       ddisp=device number for existing window to be used for dispersion plot
	# I'm not yet using a 'names' argument that might be used for labelling.
	# Output: indices of station pairs selected in the dispersion plot, 
	# and indices of individual stations selected on the coordinate plot.
	#
        color = colors()[c(12,26,32,37,53,60,70,80,84,88,94,101,116,142,366,371,376,386,392,398,400:657)]

 	n <- nrow(coords)
 	dist.mx <- Fdist(coords)
 	ix <- matrix(1.:n, n, n)
 	#
 	jx <- #
 	t(ix)
 	if(missing(dmap)) {
 		eval(call(device))
 		dmap <- # window number for the map
 		dev.cur()
 	}
 	dev.set(dmap)
 	plim <- setplot(coords, axes = T)
 	points(coords)
 	title(main = "Geographic Coordinates")
 	if(missing(ddisp)) {
 		eval(call(device))
 		ddisp <- # window number for dispersion plot
 		dev.cur()
 	}
 	dev.set(ddisp)
 	#
 	#
 	plot(dist.mx, disp.mx, main = "Spatial Dispersion: Var(Z(x)-Z(y))")
 	# counter for points identified on dispersion plot.
 	i <- 0.
 	# counter for stations identified on coordinate plot.
 	j <- 0.
 	ss <- NULL
 	pp <- NULL
 	rx <- 1.
 	while(length(rx) > 0.) {
 		cat("Enter\n 1 to identify on coordinate plot\n", 
 			"2 to identify on dispersion plot\n", 
 			"<return> to quit\n")
 		px <- 1.
 		sx <- 1.
 		rx <- scan(n = 1.)
 		if(length(rx) > 0.) {
 			if(rx == 2.) {
 				while(length(px) > 0.) {
 					dev.set(ddisp)
cat("Identify a point on the dispersion plot (left button to select; right to stop)\n")
 					px <- identify(dist.mx, disp.mx, n = 1.,plot = F)
 					if(length(px) > 0.) {
 					pp <- c(pp, px)
 					i <- i + 1.
 					points(dist.mx[pp[i]], disp.mx[pp[i]], 
                        pch = i+1, cex = 2.,col=color[3+i])
 					dev.set(dmap)
 					points(coords[c(ix[pp[i]],
 						jx[pp[i]]), ], pch = i+1, cex = 2.,col=color[3+i])
 					lines(coords[c(ix[pp[i]], jx[pp[i]]),],col=color[3+i])
 					}
 				}
 			}
 			else {
 				while(length(sx) > 0.) {
 					dev.set(dmap)
cat("Identify a site on the coordinate plot (left button to select; right to stop)\n")
 			sx <- identify(coords[, 1.], coords[, 2.], n = 1., plot = F)
 					if(length(sx) > 0.) {
 						j <- j + 1.
 						ss <- c(ss, sx)
 						points(coords[ss[j], 1.],coords[ss[j], 2.],
 							pch = 14. + j, cex = 1.5,col=color[3+j])
 						dev.set(ddisp)
 						points(dist.mx[ss[j],  ],disp.mx[ss[j],  ],
 							pch = 14. + j, cex = 1.5, col=color[3+j])
 					}
 				}
 			}
 		}
 	}
 	selpairs <- cbind(ix[pp], jx[pp])
    dev.off(dmap)
    dev.off(ddisp)
# 	return(selpairs, ss)
    list(selpairs= selpairs, ss = ss)
 }

bgrid <- function(start, xmat, coef, xlimt, iter.limit = 10., perpc = 8.)
{
	if(length(start) != 2.)
		stop(paste(substitute(start), "is not a 2 element vector"))
	if(!is.matrix(xmat))
		stop(paste(substitute(xmat), "is not a matrix"))
	if(!is.matrix(coef))
		stop(paste(substitute(coef), "is not a matrix"))
	if(missing(xlimt))
		xlimt <- c(apply(xmat[, 1.:2.], 2., range))
	if(length(xlimt) != 4.)
		stop(paste(substitute(xlimt), "is not a 4 element vector"))
	nn <- as.integer(2000.)
	start <- as.single(start)
	xylim <- as.single(xlimt)
	datlen <- as.integer(nrow(xmat))
	datain <- matrix((xmat), datlen, ncol(xmat))
	coefdim <- as.integer(ncol(coef))
	coefin <- matrix((coef), datlen + 3., coefdim)
	grid <- matrix((0.), 2., nn)
 	fldmag <- rep(0, nn)
	ngrid <- vector("integer", 1.)
	iter.limit <- as.integer(iter.limit)
	perpc <- as.integer(perpc)
	w21 <- grid
	w22 <- grid
	w1 <- fldmag
	w2 <- fldmag
	w3 <- fldmag
	Q <- .Fortran("bgrid",
		nn = as.integer(nn),
		start=as.single(start),
		xylim=as.single(xylim),
		datain,
		datlen=as.integer(datlen),
		coefin,
		coefdim=as.integer(coefdim),
		grid = grid,
		ngrid=as.integer(ngrid),
		fldmag = as.single(fldmag),
		iter.limit=as.integer(iter.limit),
		perpc=as.integer(perpc),
		w21,
		w22,
		w1 = as.single(w1),
		w2 = as.single(w2),
		w3 = as.single(w3))
	R <- list(grid = Q$grid, ngrid = Q$ngrid, fldmag = Q$fldmag, nn = Q$nn)
	if(R$nn == 999.)
		stop("internal grid buffer overflow")
	R$grid <- R$grid[, c(1.:R$ngrid)]
	R$fldmag <- R$fldmag[c(1.:R$ngrid)]
	R$fldmag[R$grid[1.,  ] == 999.] <- NA
	R$grid[R$grid == 999.] <- NA
	return(R)
}

draw <- function(data, fs = FALSE, lwidth = c(1., 1.), lcolor = c(1., 1.), cutpts,
	limits = FALSE, optlist, pts = FALSE)
{
	if(missing(data))
		stop("no data was passed to draw")
	if(!all(names(data) == c("grid", "ngrid", "fldmag", "nn")))
		stop(paste(substitute(data), "is not a valid grid object"))
	if(!is.matrix(data$grid))
		stop(paste(substitute(data), "$grid is not a matrix", sep = "")
			)
	grid1 <- data$grid[1.,  ]
	grid2 <- data$grid[2.,  ]
	if(pts)
		points(grid1, grid2)
	if(!fs) {
		lines(grid1, grid2)
		return()
	}
	if(!missing(optlist)) {
		if(missing(lwidth) && !is.null(optlist$lwidth))
			lwidth <- optlist$lwidth
		if(missing(lcolor) && !is.null(optlist$lcolor))
			lcolor <- optlist$lcolor
		if(missing(cutpts) && !is.null(optlist$cutpts))
			cutpts <- optlist$cutpts
		if(missing(limits) && !is.null(optlist$limits))
			limits <- optlist$limits
	}
	wid.cnt <- abs(wid.range <- as.integer(lwidth[2.] - lwidth[1.]))
	col.cnt <- abs(col.range <- as.integer(lcolor[2.] - lcolor[1.]))
	linetype <- (wid.cnt == 0.) & (col.cnt == 0.)
	gridlen <- data$ngrid
	fmag <- data$fldmag
	gridmiss <- (1.:gridlen)[is.na(fmag)]
	fmag[gridmiss] <- fmag[gridmiss - 1.]
	fmag.range <- range(fmag)
	if(missing(cutpts) && (missing(optlist) || is.null(optlist$cutpts))) {
		if(!linetype)
			intcnt <- max(wid.cnt, col.cnt) + 1.
		else intcnt <- 4.
		fmag.intlen <- (fmag.range[2.] - fmag.range[1.])/intcnt
		cutpts <- seq(fmag.range[1.], fmag.range[2.], fmag.intlen)
		cutpts[1.] <- 0.
		limits <- TRUE
	}
	optlist <- list(lwidth = lwidth, lcolor = lcolor, cutpts = cutpts,
		limits = limits)
	if(!limits)
		cutpts <- c(0., cutpts, fmag.range[2.])
	fmag.breaks <- cut(fmag, cutpts)
	fmag.brk1 <- c(fmag.breaks[2.:gridlen], 0.)
	intcnt <- as.integer(length(cutpts) - 1.)
	difcnt <- as.integer(intcnt - 1.)
	m <- matrix(1.:intcnt, nrow = gridlen, ncol = intcnt, byrow = T)
	lt <- matrix(FALSE, nrow = gridlen, ncol = intcnt)
	lt[fmag.breaks == m | fmag.brk1 == m] <- TRUE
	x <- matrix(rep(grid1, intcnt), ncol = intcnt)
	y <- matrix(rep(grid2, intcnt), ncol = intcnt)
	x[lt == FALSE] <- NA
	y[lt == FALSE] <- NA
	if(!linetype) {
		if(wid.cnt > 0.)
			widths <- seq(lwidth[1.], lwidth[2.], wid.range/difcnt)
		else widths <- rep(1., intcnt)
		if(col.cnt > 0.)
			colors <- seq(lcolor[1.], lcolor[2.], col.range/difcnt)
		else colors <- rep(1., intcnt)
		for(ix in 1.:intcnt)
			lines(x[, ix], y[, ix], lty = 1., lwd = widths[ix],
				col = colors[ix])
	}
  else matlines(x, y)
	return(optlist)
}

drinteg <- function(xmat, spcoef, xlimt, ind = 2., iter.limit, ncpar = 1., fs = FALSE,
	optlist)
{
	if(!is.matrix(xmat))
		stop(paste(substitute(xmat), "is not a matrix"))
	if(!is.matrix(spcoef))
		stop(paste(substitute(spcoef), "is not a matrix"))
	if(missing(xlimt))
		xlimt <- c(apply(xmat[, 1.:2.], 2., range))
	for(x1 in 1.:ncpar) {
		start <- locator(1.)
		if(length(start) > 1.) {
			start <- c(start$x, start$y)
			if(!missing(iter.limit))
				single <- integ(start, xmat, spcoef, ind, xlimt,
					iter.limit)
			else single <- integ(start, xmat, spcoef, ind, xlimt)
			if(!missing(optlist))
				optlist <- draw(single, fs = fs, optlist = 
					optlist)
			else optlist <- draw(single, fs = fs)
		}
	}
	return(single)
}


corrfit <- function(crds, Tspline, sg.fit, model = 1)
#This function estimates correlations between all the locations(new+stations)
# using the results of the SG step


#Input
#   crds : coordinates of all locations beginning with new locations
#   Tspline: the thin-spline fit from the SG-steps
#   sg.fit: the mapping resulted from the SG method
#   Model: variogram model; 1: exponential  2: gaussian
#Output
#   cor: correlation matrix among the locations
{
d.loc <- seval(crds,Tspline)
h = Fdist(t(d.loc$y))
a0 = sg.fit$variogfit$a[1]
t0 = sg.fit$variogfit$t0
if (model==1) dispfit = a0 + (2-a0)*(1- exp(- t0* h))
if (model==2) dispfit = a0 + (2-a0)*(1- exp(- t0* h^2))

corfit = (2-dispfit)/2
diag(corfit) = 1
list(cor = corfit)
}


