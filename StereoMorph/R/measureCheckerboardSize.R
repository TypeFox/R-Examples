measureCheckerboardSize <- function(corner.file, nx, ruler.file=NULL, ruler.pt.size=NULL){

	# IF INPUTS ARE TEXT FILES, READ FILES INTO MATRICES
	if(is.matrix(corner.file)){
		corners <- corner.file
	}else{
		corners <- as.matrix(read.table(corner.file))
	}
	if(!is.null(ruler.file)){
		if(is.matrix(ruler.file)){
			ruler.pts <- ruler.file
		}else{
			if(grepl('.txt$', ruler.file)){ruler.pts <- as.matrix(read.table(ruler.file, row.names=1))}else{ruler.pts <- ruler.file}
		}
	}
	
	# TRY SOME CONTROL CASES
	#corners <- 20*scale(cbind(rep(1:5, 5), c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5), rep(5, 5))))
	#rotm <- matrix(c(cos(0.1), sin(0.1), -sin(0.1), cos(0.1)), nrow=2)
	#corners <- corners %*% rotm
	#corners <- corners + rnorm(length(corners), sd=0.0)
	#plot(corners)
	#print(distancePointToPoint(corners[1, ], corners[2, ]))

	# GET OTHER INTERNAL CORNER DIMENSION FROM CORNER DIMENSIONS
	ny <- nrow(corners) / nx

	# FIT MODEL TO POINTS
	resample_grid_pts <- resampleGridImagePoints(corners, nx = nx, rx = 2, ry = 2)

	# FIND DIFFERENCE IN EDGE LENGTH
	side_lengths <- c(
		distancePointToPoint(resample_grid_pts$pts[1, ], resample_grid_pts$pts[2, ]),
		distancePointToPoint(resample_grid_pts$pts[2, ], resample_grid_pts$pts[4, ]),
		distancePointToPoint(resample_grid_pts$pts[4, ], resample_grid_pts$pts[3, ]),
		distancePointToPoint(resample_grid_pts$pts[3, ], resample_grid_pts$pts[1, ]))

	# FIT MINIMAL PARAMETER CHECKERBOARD MODEL TO POINTS
	nlminb_fit_x <- nlminb(
		start=c(corners[1, 1], corners[2, 1]-corners[1, 1], 0), 
		objective=gridPointsFitError, points=corners[, 1], nx=nx, ny=ny)
	nlminb_fit_y <- nlminb(
		start=c(corners[1, 2], corners[2, 2]-corners[1, 2], 0), 
		objective=gridPointsFitError, points=corners[, 2], nx=nx, ny=ny)

	# GET POINT FIT
	grid_pts_fit <- cbind(gridPointsFit(nlminb_fit_x$par, nx=nx, ny=ny), gridPointsFit(nlminb_fit_y$par, nx=nx, ny=ny))
	dist_grid_fit <- sqrt(rowSums((corners - grid_pts_fit)^2))

	# GET SQUARE SIZE FROM MODEL FIT
	square_size_px <- sqrt(nlminb_fit_x$par[2]^2 + nlminb_fit_y$par[2]^2)

	# PLOT POINT FIT
	#plot(corners)
	#points(grid_pts_fit, col='red', cex=0.5)

	if(is.null(ruler.file)){
		l <- list(
			side.lengths=side_lengths,
			dist.corner.fit.mean=mean(dist_grid_fit),
			dist.corner.fit.sd=sd(dist_grid_fit),
			square.size.px=square_size_px, 
			dist.ruler.fit.mean=NULL,
			dist.ruler.fit.sd=NULL,
			ruler.size.px=NULL, 
			rwu.per.px=NULL, 
			unit=NULL
			)
		class(l) <- 'measureCheckerboardSize'
		return(l)
	}

	# TRY SOME CONTROL CASES
	#ruler.pts <- cbind((1:500) + rnorm(500, sd=0.2), (1:500) + rnorm(500, sd=0.2))
	#plot(ruler.pts)
	#ruler.pts <- cbind(1:25, rep(0, 25))
	#ruler.pts <- cbind(rep(0, 25), 1:25)
	#print(ruler.pts)

	# FIT MODEL TO POINTS
	nlminb_fit_x <- nlminb(start=c(ruler.pts[1, 1], ruler.pts[2, 1]-ruler.pts[1, 1]), objective=gridPointsFitError, nx=nrow(ruler.pts), points=ruler.pts[, 1])
	nlminb_fit_y <- nlminb(start=c(ruler.pts[1, 2], ruler.pts[2, 2]-ruler.pts[1, 2]), objective=gridPointsFitError, nx=nrow(ruler.pts), points=ruler.pts[, 2])
	#print(nlminb_fit_x);print(nlminb_fit_y)

	# GET POINT FIT
	ruler_pts_fit <- cbind(gridPointsFit(nlminb_fit_x$par, nx=nrow(ruler.pts)), gridPointsFit(nlminb_fit_y$par, nx=nrow(ruler.pts)))
	dist_ruler_fit <- sqrt(rowSums(ruler.pts - ruler_pts_fit)^2)

	# PLOT POINT FIT
	#plot(ruler.pts, asp=1)
	#points(ruler_pts_fit, col='red', cex=2)

	# GET DISTANCE BETWEEN RULER POINTS IN PIXELS
	ruler_size_px <- sqrt(nlminb_fit_x$par[2]^2 + nlminb_fit_y$par[2]^2)
	#print(ruler_size_px)

	# SIMPLE INTERPOINT DISTANCE WITH RE-SAMPLING WILL OVERESTIMATES INTERPOINT DISTANCE
	#ruler_size_px <- mean(sqrt(rowSums((ruler.pts[1:(nrow(ruler.pts)-1), ] - ruler.pts[2:nrow(ruler.pts), ])^2)))
	#print(ruler_size_px)

	# IF RULER.PT.SIZE HAS UNITS, REMOVE AND SAVE UNITS
	ruler_pt_size <- as.numeric(gsub(pattern='([0-9|.])([ ]*)([A-Za-z]*)', replacement='\\1', x=ruler.pt.size))
	unit <- gsub(pattern='([0-9|.])([ ]*)([A-Za-z]*)', replacement='\\3', x=ruler.pt.size)

	# GET SQUARE SIZE IN REAL-WORLD UNITS
	square_size_rwu <- square_size_px * (ruler_pt_size / ruler_size_px)

	# GET REAL-WORLD UNITS PER PIXEL
	rwu_per_px <- square_size_rwu / square_size_px

	l <- list(
		side.lengths=side_lengths,
		dist.corner.fit.mean=mean(dist_grid_fit),
		dist.corner.fit.sd=sd(dist_grid_fit),
		square.size.px=square_size_px, 
		square.size.rwu=square_size_rwu, 
		dist.ruler.fit.mean=mean(dist_ruler_fit),
		dist.ruler.fit.sd=sd(dist_ruler_fit),
		ruler.size.px=ruler_size_px, 
		rwu.per.px=rwu_per_px, 
		unit=unit
		)
	class(l) <- 'measureCheckerboardSize'
	return(l)
}

summary.measureCheckerboardSize <- function(object, ...){

	side.lengths <- object$side.lengths

	r <- ''
	r <- c(r, '\nmeasureCheckerboardSize Summary\n')
	r <- c(r, '\tChecks for planarity of checkerboard\n')
	r <- c(r, "\t\tMean length of opposing sides: ", 
		   round(mean(c(side.lengths[1], side.lengths[3])), 2), " px x ",
		   round(mean(c(side.lengths[2], side.lengths[4])), 2), " px",
		   "\n")
	r <- c(r, "\t\tDifferences in lengths of opposing sides: ", 
		   round(side.lengths[1]-side.lengths[3], 2), " px x ",
		   round(side.lengths[2]-side.lengths[4], 2), " px",
		   "\n")

	r <- c(r, '\tSimple checkerboard model fit\n')
	r <- c(r, "\t\tSquare size in pixels: ", format(object$square.size.px), " px\n")
	r <- c(r, "\t\tMean distance of points from model: ", 
		   round(object$dist.corner.fit.mean, 2), " px +/- ", round(object$dist.corner.fit.sd, 2), "",
		   "\n")

	if(is.null(object$dist.ruler.fit.mean)){
		class(r) <- "summary.measureCheckerboardSize"
		return(r)
	}

	r <- c(r, '\tRuler points model fit\n')
	r <- c(r, "\t\tDistance between points on ruler: ", format(object$ruler.size.px), " px\n")
	r <- c(r, "\t\tMean distance of points from model: ", 
		   round(object$dist.ruler.fit.mean, 2), " px +/- ", round(object$dist.ruler.fit.sd, 2), "",
		   "\n")

	r <- c(r, '\tReal-world units\n')
	r <- c(r, "\t\tSquare size in real-world units: ", format(object$square.size.rwu), " ", object$unit,"\n")
	r <- c(r, "\t\tReal-world units per pixel: ", format(object$rwu.per.px), " ", object$unit,"\n")
	r <- c(r, '\n')

	class(r) <- "summary.measureCheckerboardSize"
	return(r)
}

print.summary.measureCheckerboardSize <- function(x, ...) cat(x, sep='')