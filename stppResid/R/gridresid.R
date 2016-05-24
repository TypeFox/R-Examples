gridresid <- function(X, cifunction, theta = NULL, lambda = NULL, grid = c(10, 10), gf = NULL, resid = c("raw", "pearson", "inverse"), algthm = c("cubature", "mc", "miser", "none"), n = 100, n.miser = 10000, tol = 1e-05, maxEval = 0, absError = 0, ints = NULL)
{
	if(!is.stpp(X))
		stop("X must be an object of type stpp")
	if(!is.null(lambda)) {
		if(length(lambda) != length(X$x))
			stop("lambda must be same length as number of points")
		lamb1 <- lambda
	} else { 
		if(is.null(theta)) {
			lamb1 <- cifunction(X)
		} else 
			lamb1 <- cifunction(X, theta)
	}	
	if(is.null(gf)) {
		gf <- make.grid(stwin(X$xcoord, X$ycoord, X$tcoord), grid)
	} else
		if(!is.stgrid(gf))
			stop("gf must be an object of type stgrid")
	if(missing(resid))
		resid = "raw"
	if(missing(algthm))
		algthm = "cubature"
	if(!is.null(ints))
		algthm = "none"
	if((algthm == "none") & is.null(ints))
		algthm = "cubature"
	countf.s <- function(xy)
	{
		which((gf[[1]]$xmin <= xy[1]) & (gf[[1]]$xmax > xy[1]) & (gf[[1]]$ymin <= xy[2]) & (gf[[1]]$ymax > xy[2]))
	}
	if(resid == "raw") {
		place <- apply(cbind(X$x, X$y), 1, countf.s)
		count <- tabulate(place)
		if(length(count) != nrow(gf$grid.full))
			count <- c(count, rep(0, nrow(gf$grid.full) - length(count)))
    if(algthm == "cubature") {
      bins <- bin.info.cubature(X, cifunction, theta = theta, gf = gf, type = 1, tol = tol, maxEval = maxEval, absError = absError)
      ints <- bins$integral
      residuals <- count - ints
      y <- list(X = X, grid = gf, residuals = residuals)
      y <- c(y, bins)
    }
		if(algthm == "mc") {
			bins <- bin.info(X, cifunction, theta = theta, gf = gf, type = 1, n.def = n)
			ints <- bins$integral
			residuals <- count - ints
			y <- list(X = X, grid = gf, residuals = residuals)
			y <- c(y, bins)
		}
		if(algthm == "miser") {
			bins <- bin.info.MISER(X, cifunction, theta = theta, gf = gf, type = 1, n = n.miser)
			ints <- bins$integral
			residuals <- count - ints
			y <- list(X = X, grid = gf, residuals = residuals)
			y <- c(y, bins)
		}
		if(algthm == "none") {
			residuals <- count - ints
			y <- list(X = X, grid = gf, residuals = residuals)
		}
	}
	if(resid == "pearson") {
		count <- rep(0, nrow(gf$grid.full))
		place <- apply(cbind(X$x, X$y), 1, countf.s)
		sums <- aggregate(cbind(place, sqrt(lamb1)), by = list(place), sum)
		count[sums[,1]] <- sums[,3]
		if(algthm == "cubature") {
		  bins <- bin.info.cubature(X, cifunction, theta = theta, gf = gf, type = 2, tol = tol, maxEval = maxEval, absError = absError)
		  ints <- bins$integral
		  residuals <- count - ints
		  y <- list(X = X, grid = gf, residuals = residuals)
		  y <- c(y, bins)
		}
		if(algthm == "mc") {
			bins <- bin.info(X, cifunction, theta = theta, gf = gf, type = 2, n.def = n)
			ints <- bins$integral
			residuals <- count - ints
			y <- list(X = X, grid = gf, residuals = residuals)
			y <- c(y, bins)
		}
		if(algthm == "miser") {
			bins <- bin.info.MISER(X, cifunction, theta = theta, gf = gf, type = 2, n = n.miser)
			ints <- bins$integral
			residuals <- count - ints
			y <- list(X = X, grid = gf, residuals = residuals)
			y <- c(y, bins)
		}
		if(algthm == "none") {
			residuals <- count - ints
			y <- list(X = X, grid = gf, residuals = residuals)
		}
	}
	if(resid == "inverse") {
		count <- rep(0, nrow(gf$grid.full))
		vol <- diff(X$xcoord) / grid[1] * diff(X$ycoord) / grid[2] * diff(X$tcoord) 
		place <- apply(cbind(X$x, X$y), 1, countf.s)
		sums <- aggregate(cbind(place, lamb1^(-1)), by = list(place), sum)
		count[sums[,1]] <- sums[,3]
		residuals <- count - rep(vol, length(count))
		y <- list(X = X, grid = gf, residuals = residuals)
	}
	class(y) <- "gridresid"
	return(y)
}