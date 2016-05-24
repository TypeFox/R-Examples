devresid <- function(X, cifunction1, cifunction2, theta1 = NULL, theta2 = NULL, lambda1 = NULL, lambda2 = NULL, grid = c(10, 10), gf = NULL, algthm1 = c("cubature", "mc", "miser", "none"), algthm2 = c("cubature", "mc", "miser", "none"), n = 100, n1.miser = 10000, n2.miser = 10000, tol = 1e-05, maxEval = 0, absError = 0, ints1 = NULL, ints2 = NULL)
{
	if(!is.stpp(X))
		stop("X must be an object of type stpp")
	if(!is.null(lambda1)) {
		if(length(lambda1) != length(X$x))
			stop("lambda1 must be same length as number of points") 
		lamb1 <- lambda1
	} else { 
		if(is.null(theta1)) {
			lamb1 <- cifunction1(X)
		} else 
			lamb1 <- cifunction1(X, theta1)
	}
	if(!is.null(lambda2)) { 
		if(length(lambda2) != length(X$x))
			stop("lambda2 must be same length as number of points")
		lamb2 <- lambda2
	} else { 
		if(is.null(theta2)) {
			lamb2 <- cifunction2(X)
		} else 
			lamb2 <- cifunction2(X, theta2)
	}	
	if(is.null(gf)) {
		gf <- make.grid(stwin(X$xcoord, X$ycoord, X$tcoord), grid)
	} else
		if(!is.stgrid(gf))
			stop("gf must be an object of type stgrid")
  	if(missing(algthm1))
  		algthm1 = "cubature"
  	if(missing(algthm2))
  		algthm2 = "cubature"
  	if(!is.null(ints1)) 
		algthm1 = "none"
	if(!is.null(ints2))
		algthm2 = "none"
	if((algthm1 == "none") & (is.null(ints1)))
		algthm1 = "cubature"
	if((algthm2 == "none") & (is.null(ints2)))
		algthm2 = "cubature"
	countf.s <- function(xy)
	{
		which((gf[[1]]$xmin <= xy[1]) & (gf[[1]]$xmax > xy[1]) & (gf[[1]]$ymin <= xy[2]) & (gf[[1]]$ymax > xy[2]))
	}
	count1 <- rep(0, nrow(gf$grid.full))
	count2 <- rep(0, nrow(gf$grid.full))
	place <- apply(cbind(X$x, X$y), 1, countf.s)
	sums <- aggregate(cbind(place, log(lamb1), log(lamb2)), by = list(place), sum)
	count1[sums[,1]] <- sums[,3]
	count2[sums[,1]] <- sums[,4]
  if(algthm1 == "cubature") {
    cat("Model 1 integrals (cubature):\n")
    bins1 <- bin.info.cubature(X, cifunction = cifunction1, theta = theta1, gf = gf, type = 1, tol = tol, maxEval = maxEval, absError = absError)
    int1 <- bins1$integral
    names(bins1) <- c("n.1", "integral.1", "sd.lambda.1")
  }
	if(algthm2 == "cubature") {
	  cat("Model 2 integrals (cubature):\n")
	  bins2 <- bin.info.cubature(X, cifunction = cifunction2, theta = theta2, gf = gf, type = 1, tol = tol, maxEval = maxEval, absError = absError)
	  int2 <- bins2$integral
	  names(bins2) <- c("n.2", "integral.2", "sd.lambda.2")
	}
  	if(algthm1 == "miser") {
  		cat("Model 1 integrals (miser with ", n1.miser, " points):\n")
    	bins1 <- devbin.info.MISER(X, cifunction1, theta = theta1, gf = gf, n = n1.miser)
    	int1 <- bins1$integral
  		names(bins1) <- c("n.1", "integral.1", "mean.lambda.1", "sd.lambda.1", "app.pts.1")
  	}
  	if(algthm2 == "miser") {
  		cat("Model 2 integrals (miser with ", n2.miser, " points):\n")
  		bins2 <- devbin.info.MISER(X, cifunction2, theta = theta2, gf = gf, n = n2.miser)
  		int2 <- bins2$integral
  		names(bins2) <- c("n.2", "integral.2", "mean.lambda.2", "sd.lambda.2", "app.pts.2")
  	}
  	if(algthm1 == "mc") {
  		cat("Model 1 integrals (mc):\n")
  		bins1 <- bin.info(X, cifunction = cifunction1, theta = theta1, gf = gf, type = 1, n.def = n)
  		int1 <- bins1$integral
  		names(bins1) <- c("n.1", "integral.1", "mean.lambda.1", "sd.lambda.1")
  	}
  	if(algthm2 == "mc") {
  		cat("Model 2 integrals (mc):\n")
  		bins2 <- bin.info(X, cifunction = cifunction2, theta = theta2, gf = gf, type = 1, n.def = n)
  		int2 <- bins2$integral
  		names(bins2) <- c("n.2", "integral.2", "mean.lambda.2", "sd.lambda.2")
  	}
  	if(algthm1 == "none") {
  		int1 <- ints1
  		bins1 <- c()
  	}
  	if(algthm2 == "none") {
  		int2 <- ints2
  		bins2 <- c()
  	}
    residuals <- (count1 - int1) - (count2 - int2)
  	y <- list(X = X, grid = gf, residuals = residuals)
  	y <- c(y, bins1, bins2)
  	class(y) <- "devresid"
  	return(y)
}