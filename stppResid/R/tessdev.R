tessdev <- function(X, cifunction1, cifunction2, theta1 = NULL, theta2 = NULL, lambda1 = NULL, lambda2 = NULL, algthm1 = c("mc", "miser", "none"), algthm2 = c("mc", "miser", "none"), n = 100, n1.miser = 10000, n2.miser = 10000, ints1 = NULL, ints2 = NULL)
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
	data <- data.frame(cbind(X$x, X$y))
	if(missing(algthm1))
  		algthm1 = "miser"
  	if(missing(algthm2))
  		algthm2 = "miser"
  	if(!is.null(ints1)) 
		algthm1 = "none"
	if(!is.null(ints2))
		algthm2 = "none"
	if((algthm1 == "none") & (is.null(ints1)))
		algthm1 = "miser"
	if((algthm2 == "none") & (is.null(ints2)))
		algthm2 = "miser"
	if(nrow(unique(data)) != length(X$x)) {
		cat("Warning message: \nOverlapping points\n")
		Y <- unique(data)
		x.tess <- Y[,1]
		y.tess <- Y[,2]
		rw <- c(X$xcoord, X$ycoord)
		vor <- deldir(x.tess, y.tess, rw = rw, digits=20)
		areas <- vor$summary[,8]
		tl <- tile.list(vor)
		if(algthm1 == "miser") {
			cat("Model 1 integrals (miser with ", n1.miser, " points):\n")
			tess1 <- tess.info.MISER(X, cifunction1, theta = theta1, areas, tl, n1.miser)
			int1 <- tess1$integral
		}
		if(algthm2 == "miser") {
			cat("Model 2 integrals (miser with ", n2.miser, "points):\n")
			tess2 <- tess.info.MISER(X, cifunction2, theta = theta2, areas, tl, n2.miser)
			int2 <- tess2$integral
		}
		if(algthm1 == "none") {
			int1 <- ints1
			tess1 <- c()
		}
		if(algthm2 == "none") {
			int2 <- ints2
			tess2 <- c()
		}
		if(algthm1 == "mc") {
			cat("Model 1 integrals (mc):\n")
			tess1 <- tess.info(X, cifunction1, theta = theta1, areas, tl, n.def = n)
			int1 <- tess1$integral
		}
		if(algthm2 == "mc") {
			cat("Model 2 integrals (mc):\n")
			tess2 <- tess.info(X, cifunction2, theta = theta2, areas, tl, n.def = n)
			int2 <- tess2$integral
		}
		sum.l <- function(tl) {
			cells <- as.points(tl$x, tl$y)
			place <- inpip(data, cells)
			total1 <- sum(lamb1[place])
			total2 <- sum(lamb2[place])
			return(c(total1, total2))
		}
		tot.l <- sapply(tl, sum.l)
		residuals <- (tot.l[1,] - int1) - (tot.l[2,] - int2)
	} else {
		x.tess <- X$x
		y.tess <- X$y
		rw <- c(X$xcoord, X$ycoord)
		vor <- deldir(x.tess, y.tess, rw = rw, digits=20)
		areas <- vor$summary[,8]
		tl <- tile.list(vor)
		if(algthm1 == "miser") {
			cat("Model 1 integrals (miser with ", n1.miser, " points):\n")
			tess1 <- tess.info.MISER(X, cifunction1, theta = theta1, areas, tl, n1.miser)
			int1 <- tess1$integral
		}
		if(algthm2 == "miser") {
			cat("Model 2 integrals (miser with ", n2.miser, "points):\n")
			tess2 <- tess.info.MISER(X, cifunction2, theta = theta2, areas, tl, n2.miser)
			int2 <- tess2$integral
		}
		if(algthm1 == "none") {
			int1 <- ints1
			tess1 <- c()
		}
		if(algthm2 == "none") {
			int2 <- ints2
			tess2 <- c()
		}
		if(algthm1 == "mc") {
			cat("Model 1 integrals (mc):\n")
			tess1 <- tess.info(X, cifunction1, theta = theta1, areas, tl, n.def = n)
			int1 <- tess1$integral
		}
		if(algthm2 == "mc") {
			cat("Model 2 integrals (mc):\n")
			tess2 <- tess.info(X, cifunction2, theta = theta2, areas, tl, n.def = n)
			int2 <- tess2$integral
		}
		residuals <- (lamb1 - int1) - (lamb2 - int2)
	}
	y <- list(X = X, tile.list = tl, residuals = residuals)
	y <- c(y, tess1, tess2)
	class(y) <- "tessdev"
	return(y)
}