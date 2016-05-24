# ===========================================================================
# File: "main.R"
#                        Created: 2011-01-11 12:47:48
#              Last modification: 2013-11-22 11:38:53
# Author: Bernard Desgraupes
# e-mail: <bernard.desgraupes@u-paris10.fr>
# ===========================================================================
# Plot conics given the symmetric matrix.


## 
 # ------------------------------------------------------------------------
 # 
 # "conicPlot(x, etc...)" --
 # 
 # ------------------------------------------------------------------------
 ##
conicPlot <- function(x, type='l', npoints=100, 
	sym.axes=FALSE, center=FALSE, asymptotes=FALSE, add=FALSE, 
	xlim=NULL, ylim=NULL, ax.lty=1, ax.col=palette()[1], as.lty=1, as.col=palette()[1], ...) 
{
	m <- conicGetMatrix(x)
	D <- conicCheckDet(m)
	d <- conicCheckDet(m[1:2,1:2])
	
	if (D != 0) {
		if (d > 0) {
			res <- plotEllipsis(m, D, d, npoints=npoints, 
				sym.axes=sym.axes, center=center, add=add, 
				xlim=xlim, ylim=ylim, type=type, 
				ax.lty=ax.lty, ax.col=ax.col, ...)
		} else if (d < 0) {
			res <- plotHyperbola(m, D, d, npoints=npoints, 
				sym.axes=sym.axes, center=center, asymptotes=asymptotes, add=add, 
				xlim=xlim, ylim=ylim, type=type, 
				ax.lty=ax.lty, ax.col=ax.col, as.lty=as.lty, as.col=as.col, ...)
		} else {
			res <- plotParabola(m, npoints=npoints, 
				sym.axes=sym.axes, add=add, 
				xlim=xlim, ylim=ylim, type=type, 
				ax.lty=ax.lty, ax.col=ax.col, ...)
		}
	} else {
		if (options("warn")[[1]]) {
			warning("conic is degenerate")
		}
		res <- plotDegenerateConic(m, d, npoints=npoints, 
			sym.axes=sym.axes, center=center, add=add, 
			xlim=xlim, ylim=ylim, type=type, 
			ax.lty=ax.lty, ax.col=ax.col, ...)
	}
	
	invisible(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "plotEllipsis(m)" --
 # 
 # This is the case det(B) > 0. Both eigenvalues have the same sign: this
 # must be the opposite of the sign of det(m), otherwise the conic is empty.
 # 
 # ------------------------------------------------------------------------
 ##
plotEllipsis <- function(m, D, d, npoints, 
						sym.axes, center, add, 
						xlim, ylim, type, 
						ax.lty, ax.col, ...) 
{
	# Find the eigenvalues and eigenvectors of matrix B
	B <- m[1:2,1:2]
	vp <- eigen(B)
	val <- vp$values
	vec <- vp$vectors
	r <- D/d
	res <- list(kind="ellipse")
	res[["axes"]] <- vec
	
	# D and the eignevalues must not have same sign
	empty <- (val[1] > 0 && D > 0) || (val[1] < 0 && D < 0)
	
	# Find the center
	C <- conicCenter(m)
	res[["center"]] <- C
	
	if (!empty) {
		# Precalculate the coefficients
		a <- sqrt(abs(r/val[1]))
		b <- sqrt(abs(r/val[2]))
		ct <- vec[1,1]
		st <- vec[2,1]
		
		# Build vectors of coordinates
		t <- seq(0,2*pi,len=npoints)
		x <- a*ct*cos(t)-b*st*sin(t)+C[1]
		y <- a*st*cos(t)+b*ct*sin(t)+C[2]
		
		# Plot the parametric curve
		if (type=="n") {
			res[["points"]] <- cbind(x,y)
		} else {
			if (add) {
				lines(x, y, type=type, xlim=xlim, ylim=ylim, ...)
			} else {
				plot(x, y, type=type, xlim=xlim, ylim=ylim, ...)
			}
		}
		
		
		# Draw the axes if required
		if (sym.axes && type!="n") {
			conicPlotAxes(C, ct, st, ax.lty, ax.col)
		}
		
		if (center && type!="n") {
			points(C[1],C[2])
		}
		
		# Calculate the coords of the vertices
		XS <- c(a,0,-a,0)
		YS <- c(0,b,0,-b)
		xs <- ct*XS-st*YS+C[1]
		ys <- st*XS+ct*YS+C[2]
		res[["vertices"]] <- list(x=xs, y=ys)
		
		# Calculate the coords of the foci
		df <- sqrt(abs(a^2-b^2))
		if (a>b) {
			XF <- c(df,-df)
			YF <- c(0,0)
		} else {
			XF <- c(0,0)
			YF <- c(df,-df)
		}
		xf <- ct*XF-st*YF+C[1]
		yf <- st*XF+ct*YF+C[2]
		res[["foci"]] <- list(x=xf, y=yf)
	
		# Calculate the eccentricity
		res[["eccentricity"]] <- df/max(a,b)
	
	} else {
		if (options("warn")[[1]]) {
			warning("the conic is empty")
		}
	}

	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "plotHyperbola(m)" --
 # 
 # This is the case det(B) < 0. The eigenvalues have opposite signs: the
 # position of the branches wrt the asymptotes depends on the sign of the
 # quantity -det(m)/(lambda1*det(B)).
 # 
 # ------------------------------------------------------------------------
 ##
plotHyperbola <- function(m, D, d, npoints=100, 
						sym.axes, center, asymptotes, add, 
						xlim, ylim, type, 
						ax.lty, ax.col, as.lty, as.col, ...) 
{
	# Find the eigenvalues and eigenvectors of matrix B
	B <- m[1:2,1:2]
	vp <- eigen(B)
	val <- vp$values
	vec <- vp$vectors
	r <- D/d
	res <- list(kind="hyperbola")
	res[["axes"]] <- vec

	# Precalculate the transformation matrix
	ct <- vec[1,1]
	st <- vec[2,1]

	# Find the center
	C <- conicCenter(m)
	res[["center"]] <- C
	
	# Find the center
	asy <- conicAsymptotes(m)
	res[["asymptotes"]] <- asy
	
	# Build vectors of coordinates for both branches. The extremity points,
	# corresponding to pi/2 i-e null cosinus, are removed.
	if (r/val[1] < 0) {
		a <- sqrt(-r/val[1])
		b <- sqrt(r/val[2])
	
		t <- seq(-pi/2,pi/2,len=npoints)
		t <- t[2:(npoints-1)]
		x1 <- a*ct/cos(t)-b*st*tan(t)+C[1]
		y1 <- a*st/cos(t)+b*ct*tan(t)+C[2]
		
		t <- seq(pi/2,3*pi/2,len=npoints)
		t <- t[2:(npoints-1)]
		x2 <- a*ct/cos(t)-b*st*tan(t)+C[1]
		y2 <- a*st/cos(t)+b*ct*tan(t)+C[2]

		# Vertices
		xs <- c(a*ct+C[1],-a*ct+C[1])
		ys <- c(a*st+C[2],-a*st+C[2])

		# Foci
		df <- sqrt(a^2+b^2)
		xf <- c(df*ct+C[1],-df*ct+C[1])
		yf <- c(df*st+C[2],-df*st+C[2])
		ecc <- df/a

		# Suggested plotting ranges
		ar <- 3*a
		Rx <- c(ar*ct+C[1],-ar*ct+C[1])
		Ry <- c(ar*st+C[2],-ar*st+C[2])
	} else {
		a <- sqrt(r/val[1])
		b <- sqrt(-r/val[2])
		
		t <- seq(-pi/2,pi/2,len=npoints)
		t <- t[2:(npoints-1)]
		x1 <- a*ct*tan(t)-b*st/cos(t)+C[1]
		y1 <- a*st*tan(t)+b*ct/cos(t)+C[2]

		t <- seq(pi/2,3*pi/2,len=npoints)
		t <- t[2:(npoints-1)]
		x2 <- a*ct*tan(t)-b*st/cos(t)+C[1]
		y2 <- a*st*tan(t)+b*ct/cos(t)+C[2]
		
		# Vertices
		xs <- c(-b*st+C[1],b*st+C[1])
		ys <- c(b*ct+C[2],-b*ct+C[2])
		
		# Foci
		df <- sqrt(a^2+b^2)
		xf <- c(-df*st+C[1],df*st+C[1])
		yf <- c(df*ct+C[2],-df*ct+C[2])
		ecc <- df/b

		# Suggested plotting ranges
		br <- 3*b
		Rx <- c(-br*st+C[1],br*st+C[1])
		Ry <- c(br*ct+C[2],-br*ct+C[2])
	}

	# Set ranges
	if (is.null(xlim)) {		
		xlim <- sort(Rx)
	}
	if (is.null(ylim)) {
		ylim <- sort(Ry)
	}

	# Plot the parametric curve
	if (type=="n") {
		res[["points"]] <- cbind(c(x1,x2),c(y1,y2))
	} else {
		if (add) {
			lines(x1, y1, type=type, xlim=xlim, ylim=ylim, ...)
			lines(x2, y2, type=type, xlim=xlim, ylim=ylim, ...)
		} else {
			plot(x1, y1, type=type, xlim=xlim, ylim=ylim, ...)
			lines(x2, y2, type=type, xlim=xlim, ylim=ylim, ...)
		}
	}
	
	
	if (center && type!="n") {
		points(C[1],C[2])
	}
	if (sym.axes && type!="n") {
		conicPlotAxes(C, ct, st, ax.lty, ax.col)
	}
	if (asymptotes && type!="n") {
		conicPlotAsymptotes(C, asy, as.lty, as.col)
	}
	
	# Coords of the vertices
	res[["vertices"]] <- list(x=xs, y=ys)
	
	# Coords of the foci
	res[["foci"]] <- list(x=xf, y=yf)

	# Eccentricity
	res[["eccentricity"]] <- ecc
	
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "plotParabola(m)" --
 # 
 # This is the case det(B) = 0. One of the eigenvalues is null and there is
 # no center. Coefficients a and c have the same sign.
 # 
 # ------------------------------------------------------------------------
 ##
plotParabola <- function(m, npoints=100, 
						sym.axes, add, xlim=xlim, ylim=ylim, type, 
						ax.lty, ax.col, as.lty, as.col, ...) 
{
	# Precalculate the coefficients
	a <- m[1,1]
	b <- m[1,2]
	c <- m[2,2]
	d <- m[1,3]
	e <- m[2,3]
	f <- m[3,3]
	nrm <- sqrt(a^2+b^2)
	den <- e*a-d*b
	A <- ((a+c)*nrm)/(2*den)
	B <- (d*a+e*b)/den
	C <- (f*nrm)/(2*den)
	S2 <- -B/(2*A)
	
	res <- list(kind="parabola")
	res[["axes"]] <- matrix(c(b,-a,a,b),nrow=2)/nrm
	res[["asymptotes"]] <- -b/c

	# Build vectors of coordinates in the (V_1,V_2) basis around the vertex
	rad <- sqrt(10/abs(A))
	y2 <- seq(S2-rad,S2+rad,len=npoints)
	y1 <- A*y2^2+B*y2+C

	# Convert to coordinates in the canonic basis
	x1 <- (b*y1+a*y2)/nrm
	x2 <- (-a*y1+b*y2)/nrm

	# Set the intervals depending on the endpoints
	r2 <-  range(y2)
	r1 <- A*r2^2+B*r2+C
	e1 <- (b*r1+a*r2)/nrm
	e2 <- (-a*r1+b*r2)/nrm
	
	if (is.null(xlim)) {		
		if ( max(e1) < max(x1) ) {
			xlim <- c(max(e1), max(x1))
		} else if ( min(x1) < min(e1) ) {
			xlim <- c(min(x1), min(e1))
		} else {
			xlim <- range(x1)
		}
	}
	if (is.null(ylim)) {
		if ( max(e2) < max(x2) ) {
			ylim <- c(max(e2), max(x2))
		} else if ( min(x2) < min(e2) ) {
			ylim <- c(min(x2), min(e2))
		} else {
			ylim <- range(x2)
		}
	}

	# Plot the parametric curve
	if (type=="n") {
		res[["points"]] <- cbind(x1,x2)
	} else {
		if (add) {
			lines(x1, x2, type=type, ...)
		} else {
			plot(x1, x2, type=type, xlim=xlim, ylim=ylim, ...)
		}
	}
	
	
	# Find coords of the vertex in canonic basis
	S1 <- A*S2^2+B*S2+C
	xs <- (b*S1+a*S2)/nrm
	ys <- (-a*S1+b*S2)/nrm
	res[["vertices"]] <- list(x=xs, y=ys)
	if (sym.axes && type!="n") {
		conicPlotAxes(c(xs,ys), b, -a, ax.lty, ax.col)
	}
	
	# Coords of the focus
	XF <- S2
	YF <- S1 + 1/(4*A)
	xf <- (b*YF+a*XF)/nrm
	yf <- (-a*YF+b*XF)/nrm
	res[["foci"]] <- list(x=xf, y=yf)
	
	# Eccentricity
	res[["eccentricity"]] <- 1
	
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "plotDegenerateConic(m)" --
 # 
 # This is the case where det(m)=0. In that case, the conic is a pair of
 # lines.
 # 
 # ------------------------------------------------------------------------
 ##
plotDegenerateConic <- function(m, dt, npoints=100, 
						sym.axes, center, add, 
						xlim, ylim, type, 
						ax.lty, ax.col, ...) 
{
	B <- m[1:2,1:2]
	a <- m[1,1]
	b <- m[1,2]
	c <- m[2,2]
	d <- m[1,3]
	res <- list(kind="lines")

	# Look for the points at infinity
	if (dt < 0) {
		# Find the center
		C <- conicCenter(m)
		res[["center"]] <- C
		
		if (is.null(xlim)) {
			xlim <- c(C[1]-10,C[1]+10)
		}
		if (is.null(ylim)) {
			ylim <- c(C[2]-10,C[2]+10)
		}
	
		if (add == FALSE) {
			plot(0,0, type='n', xlim=xlim, ylim=ylim,...)
		}
		
		# Draw the lines as asymptotes
		s <- conicAsymptotes(m)
		conicPlotAsymptotes(C,s)
		res[["asymptotes"]] <- s
		
		if (center) {
			points(C[1],C[2])
		}
		vp <- eigen(B)
		vec <- vp$vectors
		ct <- vec[1,1]
		st <- vec[2,1]
		res[["axes"]] <- vec
		if (sym.axes) {
			conicPlotAxes(C, ct, st, ax.lty, ax.col)
		}
	} else if (dt == 0) {
		# Slope of unique asymptotic direction
		t <- -b/c
		res[["asymptotes"]] <- -b/c
		res[["axes"]] <- c(c,-b)
		
		# The intercepts are the solutions of f(0,x_2)=0
		e <- m[2,3]
		f <- m[3,3]
		disc <- e^2-f*c
		
		if (disc > 0) {
			# Intercepts
			yi <- (-e + c(-1,1)*sqrt(disc))/c
			res[["intercepts"]] <- yi
			
			if (is.null(xlim)) {
				xlim <- c(-10,10)
			}
			if (is.null(ylim)) {
				ylim <- c(min(0,min(yi)),max(yi)+2)
			}
			if (add == FALSE) {
				# Start an empty plot
				plot(0,0, type='n', xlim=xlim, ylim=ylim)
			}
			abline(yi[1],t)
			abline(yi[2],t)
			if (sym.axes) {
				abline(mean(yi),t, lty=ax.lty, col=ax.col)
			}
		} else {
			# This is the case of the double line
			if (c != 0) {
				res[["intercepts"]] <- -e/c
				if (is.null(xlim)) {
					xlim <- c(-2,2)
				}
				if (is.null(ylim)) {
					ylim <- c(-e/c-2,-e/c+2)
				}
				if (add == FALSE) {
					# Start an empty plot
					plot(0,0, type='n', xlim=xlim, ylim=ylim)
				}
				abline(-e/c,t)
			} else if (a != 0) {
				if (is.null(xlim)) {
					xlim <- c(-d/a-2,-d/a+2)
				}
				if (is.null(ylim)) {
					ylim <- c(-2,2)
				}
				if (add == FALSE) {
					# Start an empty plot
					plot(0,0, type='n', xlim=xlim, ylim=ylim)
				}
				abline(v=-d/a)
			}
		}
	}
	
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "conicMatrix(m)" --
 # 
 # The 'v' argument is a 6-length vector containing the
 # coefficients of a quadratic polynomial of the form:
 #     v_1 x_1^2 + v_2 x_1*x_2 + v_3 x_2^2 + v_4 x_1 + v_5 x_2 + v_6 
 # 
 # Return the symmetric 3x3 matrix corresponding to the given vector. 
 # 
 # ------------------------------------------------------------------------
 ##
conicMatrix <- function(v) {
	v <- as.vector(v)
	if (!is.numeric(v) || length(v) != 6) {
		stop("x should be a 6-length numeric vector")
	}
	m <- matrix(nrow=3,ncol=3)
	m[1,1] <- v[1]
	m[2,2] <- v[3]
	m[3,3] <- v[6]
	m[1,2] <- m[2,1] <- v[2]/2
	m[1,3] <- m[3,1] <- v[4]/2
	m[2,3] <- m[3,2] <- v[5]/2
	return(m)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "conicCenter(x)" --
 # 
 # The 'x' argument is either a 6-length vector or a symmetric 3x3 matrix.
 # The function will raise an error if the conic has no center.
 # 
 # ------------------------------------------------------------------------
 ##
conicCenter <- function(x) {
	m <- conicGetMatrix(x)
	B <- m[1:2,1:2]
	C <- solve(B,-m[1:2,3])
	return(C)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "conicAxes(x)" --
 # 
 # The 'x' argument is either a 6-length vector or a symmetric 3x3 matrix. 
 # 
 # ------------------------------------------------------------------------
 ##
conicAxes <- function(x) {
	m <- conicGetMatrix(x)
	vp <- eigen(m[1:2,1:2])
	vec <- vp$vectors
	return(vp$vectors)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "conicAsymptotes(x)" --
 # 
 # The 'x' argument is either a 6-length vector or a symmetric 3x3 matrix. 
 # 
 # ------------------------------------------------------------------------
 ##
conicAsymptotes <- function(x) {
	m <- conicGetMatrix(x)
	c <- m[2,2]
	b <- m[1,2]
	a <- m[1,1]
	res <- vector(mode="numeric")
	dt <- b^2-a*c
	if (dt >= 0) {
		if (c != 0) {
			res <- (-b+c(-1,1)*sqrt(dt))/c
		} else {
			res <- c(-a/(2*b),Inf)
		}
	}
	
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "conicPlotAxes(m)" --
 # 
 # C is the center, ct and st are respectively the cosinus and the sinus of
 # the angle between the vectors V_1 and e_1.
 # 
 # ------------------------------------------------------------------------
 ##
conicPlotAxes <- function(C, ct, st, ax.lty=1, ax.col=palette()[1]) {
	if (ct == 0 || st == 0) {
		abline(h=C[2], col=ax.col, lty=ax.lty)
		abline(v=C[1], col=ax.col, lty=ax.lty)
	} else {
		tg <- st/ct
		abline(C[2]-tg*C[1], tg, col=ax.col, lty=ax.lty)
		ctg <- ct/st
		abline(C[2]+ctg*C[1], -ctg, col=ax.col, lty=ax.lty)
	}
}


## 
 # ------------------------------------------------------------------------
 # 
 # "conicPlotAsymptotes(C, s)" --
 # 
 # C is the center, s is 2-length vector containing the slopes of the
 # asymptotes.
 # 
 # ------------------------------------------------------------------------
 ##
conicPlotAsymptotes <- function(C, s, as.lty=1, as.col=palette()[1]) {
	for (i in 1:2) {
		abline(C[2]-s[i]*C[1],s[i], col=as.col, lty=as.lty)
	}
}


## 
 # ------------------------------------------------------------------------
 # 
 # "conicGetMatrix(x)" --
 # 
 # x can be either a 6-length vector or a 3x3 symmetric matrix. The
 # function checks the validity and returns the symmetric matrix
 # representing the conic.
 # 
 # ------------------------------------------------------------------------
 ##
conicGetMatrix <- function(x) {
	if (is.vector(x)) {
		m <- conicMatrix(x)
	} else if (is.matrix(x)) {
		d <- dim(x)
		if (!is.numeric(x)) {
			stop("x should be a numeric matrix")
		}
		if (d[1] != 3 || d[2] != 3) {
			stop("x should be a 3x3 matrix")
		}
		if (!all(x==t(x))) {
			stop("x should be a symmetric matrix")
		}
		m <- x
	} else {
		stop("x should be a 6-length vector or a 3x3 symmetric matrix")	
	}	

	return(m)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "conicCheckDet(x)" --
 # 
 # Compute a determinant. If it is less than 1e-14, replace it by 0 and
 # emit a warning.
 # 
 # ------------------------------------------------------------------------
 ##
conicCheckDet <- function(m, emit=TRUE) {
	dt <- det(m)
	if (abs(dt)<1e-14) {
		if (emit && options("warn")[[1]]) {
			warning("determinant less than 1e-14, replaced by 0. Results may be inaccurate.")
		}
		dt <- 0
	}

	return(dt)
}
