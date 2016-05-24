# calculate resistance for a truncated conical conductor
# dimensions in microns
# resistivity in ohm.cm
# default rho is for 150 mM KCl
# result in ohms
resist.cone <- function ( l, r1, r2, rho=51 )
{
	# convert microns to cm
	l <- l * 1e-4;
	r1 <- r1 * 1e-4;
	r2 <- r2 * 1e-4;
	
	return ( rho * l * abs((1/r1) - (1/r2)) / (pi * abs(r1-r2)) );
}

# access resistance to hole of radius r
# calculated according to Hall 1975
# r in microns
# resistivity in ohm.cm
# result in ohms
resist.access <- function ( r, rho=51 )
{
	return ( rho / (4e-4 * r) )
}

# combined series and access resistances for a single conical
# segment of given length, half-cone angle and tip radius
# (distances in microns, angle in radians)
resist.total <- function ( r, l, theta, rho=51 )
{
	r2 <- r + l * tan(theta)
	return (resist.cone(l, r, r2, rho) + resist.access(r, rho))
}

# breakage resistance based on the above, for fitting
# result is converted results to Mohm
resist.breakage <- function ( x, theta, r, rho=51, l=1000 )
{
	result <- NA
	if ( theta > 0 & theta < pi/8 & r > 0 & r < l )
	{
		result <- resist.total(r + x * tan(theta),l-x,theta,rho) * 1e-6
	}
	return(result)
}

# residual sum of squares error function based on the above, for minimization
err.breakage <- function ( params, data, rho=51, l=1000 )
{
	y.pred <- resist.breakage( data$x, params[1], params[2], rho, l )
	y.actual <- data$y
	y.diff <- y.pred - y.actual
	SS <- sum(y.diff * y.diff)
	return(SS)
}

# fit estimated resistance-breakage data
# data is an xy list, where x is distance
# (or else a table including Mohm and Z)
# will plot result if do.plot is TRUE (default)
fit.breakage <- function ( data, start=list(theta=3*pi/180, r=0.05), rho=51, l=1000, do.plot=TRUE, ... )
{
	# locally rename columns if data is supplied as Z and Mohm
	if ( length(names(data)[names(data)=="x"]) == 0 )
	{
		names(data)[names(data)=="Z"] <- "x"
		names(data)[names(data)=="Mohm"] <- "y"
	}
	
	opt <- optim(start,
				 function(p) {err.breakage(p, data, rho, l)},
				 lower=c(0.1 * pi/180, 0.001),
				 upper=c(20 * pi/180, 50),
				 hessian=TRUE,
				 method="L-BFGS-B")
	if ( opt$convergence != 1 )
	{
		theta <- opt$par[1]
		r <- opt$par[2]
		err <- opt$value
		
		cvm <- solve(opt$hessian)
		theta.sd <- sqrt(cvm[1,1])
		r.sd <- sqrt(cvm[2,2])
		
		if ( do.plot )
		{
			xx <- 0.01 * (1:(100 * max(10, ceiling(max(data$x)))))
			plot ( xx,
			       resist.breakage(xx, theta, r, rho=rho, l=l),
			       type="l",
			       ylab="Resistance (Mohm)",
			       xlab="Breakage distance (um)",
			       ylim=c(0, 1.1 * max(data$y)),
			       ... )
			points ( data, ... )
		}
		
		return ( list(theta=theta, r=r, degrees=180*theta/pi,
					  theta.sd=theta.sd, r.sd=r.sd, degrees.sd=180*theta.sd/pi,
					  err=err, opt=opt) )
	}
	else
	{
		return ( list(theta=NA, r=NA, degrees=NA, err=NA, opt=opt) )
	}
}

# plot the error function surface over a range around the fit results
# as an illustration of fit sensitivity
fit.sensitivity.plot <- function ( data,
							       fit,
							       rho=51,
							       l=1000,
							       r.range=0.015,
							       theta.range=pi/360,
							       steps=100,
							       nlevels=200,
							       r.squared=TRUE,
							       bound.at=0.99,
							       ... )
{
	# locally rename columns if data is supplied as Z and Mohm
	if ( length(names(data)[names(data)=="x"]) == 0 )
	{
		names(data)[names(data)=="Z"] <- "x"
		names(data)[names(data)=="Mohm"] <- "y"
	}
	
	z <- mat.or.vec(nr=steps, nc=steps)
	rr <- fit$r - r.range + ((1:steps) * r.range * 2 / steps)
	tt <- fit$theta - theta.range + ((1:steps) * theta.range * 2 / steps)
	
	for ( ir in 1:steps )
	{
		for ( it in 1:steps )
		{
			z[ir, it] <- err.breakage(c(tt[it], rr[ir]), data, rho, l)
		}
	}

	rownames(z) <- rr
	colnames(z) <- tt
	
	f.label <- "Error function sensitivity"
	err.label <- signif(fit$err,4)
	
	if (r.squared)
	{
		err <- data$y - mean(data$y)
		ss.err <- sum(err * err)
		
		z <- 1 - z/ss.err
		
		f.label <- "Coefficient of Determination (R^2)"
		err.label <- signif(1-fit$err/ss.err,4)
	}
	
	result <- list(x=tt, y=rr, z=z)
	
	contour(x=tt * 180/pi, y=rr * 1e3, z=z, nlevels=nlevels, xlab="Cone Angle (deg)", ylab="Tip Radius (nm)", main=f.label, ...)
	contour(x=tt * 180/pi, y=rr * 1e3, z=z, levels=bound.at, add=TRUE, col=2)
	points(x=fit$degrees, y=fit$r * 1e3, pch=16, col=2)
	text(x=fit$degrees, y=fit$r * 1e3, col=2, pos=3, labels=err.label, cex=0.6)
	
	invisible( result )
}


# pass a vector of heights
# returns a logical vector of those points at
# the bottom of a fall
# window specifies the minimum number of previous
# points descending
find.bottom <- function(x, window=50, box.size=9, clip.ends=TRUE )
{
	if ( box.size > 1 )
	{
		x <- filter(x, method="convolution", filter=rep(1,box.size)/box.size)
		x[is.na(x)] <- 0
	}
	
	d <- c(0, diff(x)) < 0
	result <- logical(length(x))
	zero.count <- 0
	for ( ii in 1:length(x) )
	{
		if ( d[ii] )
		{
			zero.count <- zero.count + 1
		}
		else
		{
			if ( zero.count >= window )
			{
				result[ii] <- TRUE
			}
			
			zero.count <- 0
		}
	}
	
	if ( clip.ends )
	{
		result[1:min(window,length(result))] <- FALSE
		result[max(1,length(result)-window):length(result)] <- FALSE
	}
	
	return(result)
}

# for each stretch of data vector x that lies *between* bottom
# points btm (as found using find.bottom above), apply function f
# and return a vector of the results
# note that unbounded stretches at each end are dropped, so there
# will be one less result than the number of bottom points
# at least one of btm and idx must be supplied
apply.breaks <- function (x, btm=NULL, f=median, idx=(1:length(btm))[btm])
{
	result <- NULL
	for ( ii in 1:(length(idx)-1) )
	{
		lo <- idx[ii] + 1
		hi <- idx[ii+1] - 1
		if ( lo <= hi )
		result <- c(result, f(x[lo:hi]))
	}
	
	return(result)
}

# plot a current-distance or resistance-distance relation for pipette breakage data
# where distance is the bottom point extension Z
# current is the median interval current pA
# resistance is calculated if voltage is supplied in the mV argument
# optionally plot a regression line
# returns a data frame of the point data
breakage.plot <- function ( x,
							time.limits = NULL,
							btm = find.bottom(x$Z),
							pch=16,
							col=1,
							lty="solid",
							col.l=2,
							plot.line=FALSE,
							f=median,
							mV=NA,
							... )
{
	idx <- (1:length(btm))[btm]
	
	if ( length(time.limits) == 2 )
	{
		idx <- idx[(x$s[idx] > time.limits[1]) & (x$s[idx] < time.limits[2])]
	}
	
	pA <- apply.breaks(x$pA, f=f, idx=idx)
	Z <- x$Z[idx[2:length(idx)]]
	
	result <- NULL
	
	if ( ! is.na(mV) )
	{
		Mohm <- mV/pA * 1e3
		
		plot(Mohm ~ Z, pch=pch, col=col, xlab="Surface location (um)", ylab="Resistance (Mohm)", ...)
	
		if ( plot.line )
		{
			abline(lm(Mohm ~ Z), lty=lty, col=col.l, ...)
		}
		
		result <- data.frame(Z=Z, pA=pA, Mohm=Mohm)
	}
	else
	{
		plot(pA ~ Z, pch=pch, col=col, xlab="Surface location (um)", ylab="I_ref (pA)", ...)
	
		if ( plot.line )
		{
			abline(lm(pA ~ Z), lty=lty, col=col.l, ...)
		}
		
		result <- data.frame(Z=Z, pA=pA)
	}
	
	invisible(result)
}

# pick point clusters from a table including Mohm and Z values
# (eg as returned by breakage.plot)
# and return a list of the median of each cluster
# uses the Imap function select.pts for the selection
# choose an empty polygon to terminate
break.clust <- function ( data, zero.invert=TRUE )
{
	plot(Mohm ~ Z, data=data)
	pts <- data[T,c("Z", "Mohm")]
	
	result <- NULL
	
	cont <- TRUE
	while ( cont )
	{
		chosen <- try(select.pts(pts), silent=TRUE)
		if ( class(chosen) == "try-error" )
		{
			cont <- FALSE
		}
		else
		{
			result <- rbind(result, data.frame(Z=median(chosen$Z), Mohm=median(chosen$Mohm)))
		}
	}
	
	if ( zero.invert )
	{
		result$Z <- max(result$Z) - result$Z
	}
		
	return(result)
}
