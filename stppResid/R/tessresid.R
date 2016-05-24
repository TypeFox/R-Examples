tessresid <- function(X, cifunction, theta = NULL, algthm = c("mc", "miser", "none"), n = 100, n.miser = 10000, ints = NULL)
{
	if(!is.stpp(X))
		stop("X must be an object of type stpp")
	data <- data.frame(cbind(X$x, X$y))
	if(missing(algthm))
		algthm = "miser"
	if(!is.null(ints))
		algthm = "none"
	if((algthm == "none") & is.null(ints))
		algthm = "miser"
	if(nrow(unique(data)) != length(X$x)) {
		cat("Warning message: \nOverlapping points\n")
		Y <- unique(data)
		x.tess <- Y[,1]
		y.tess <- Y[,2]
		rw <- c(X$xcoord, X$ycoord)
		vor <- deldir(x.tess, y.tess, rw = rw, digits=20)
		areas <- vor$summary[,8]
		tl <- tile.list(vor)
		num <- function(tl) {
			cells <- as.points(tl$x, tl$y)
			place <- inpip(data, cells)
			n2 <- length(place)
			return(n2)
		}
		if(algthm == "mc") {
			tess <- tess.info(X, cifunction, theta = theta, areas, tl, n)
			ints <- tess$integral
			num.pts <- sapply(tl, num)
			residuals <- (num.pts - ints)/sqrt(ints)
			y <- list(X = X, tile.list = tl, residuals = residuals)
			y <- c(y, tess)
		}
		if(algthm == "miser") {
			tess <- tess.info.MISER(X, cifunction, theta = theta, areas, tl, n.miser)
			ints <- tess$integral
			num.pts <- sapply(tl, num)
			residuals <- (num.pts - ints)/sqrt(ints)
			y <- list(X = X, tile.list = tl, residuals = residuals)
			y <- c(y, tess)
		}
		if(algthm == "none") {
			num.pts <- sapply(tl, num)
			residuals <- (num.pts - ints)/sqrt(ints)
			y <- list(X = X, tile.list = tl, residuals = residuals)
		} 
	} else {
		x.tess <- X$x
		y.tess <- X$y
		rw <- c(X$xcoord, X$ycoord)
		vor <- deldir(x.tess, y.tess, rw = rw, digits=20)
		areas <- vor$summary[,8]
		tl <- tile.list(vor)
		if(algthm == "mc") {
			tess <- tess.info(X, cifunction, theta = theta, areas, tl, n)
			ints <- tess$integral
			residuals <- (1 - ints)/sqrt(ints)
			y <- list(X = X, tile.list = tl, residuals = residuals)
			y <- c(y, tess)
		}
		if(algthm == "miser") {
			tess <- tess.info.MISER(X, cifunction, theta = theta, areas, tl, n.miser)
			ints <- tess$integral
			residuals <- (1 - ints)/sqrt(ints)
			y <- list(X = X, tile.list = tl, residuals = residuals)
			y <- c(y, tess)
		}
		if(algthm == "none") {
			residuals <- (1 - ints)/sqrt(ints)
			y <- list(X = X, tile.list = tl, residuals = residuals)
		}
	}
	class(y) <- "tessresid"
	return(y)
}