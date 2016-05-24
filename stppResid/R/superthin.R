superthin<-function(X, cifunction, theta = NULL,  k = NULL, lambda = NULL)
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
	if(is.null(k)) {
		k <- length(X$x)/(diff(X$xcoord)*diff(X$ycoord)*diff(X$tcoord))
		cat("No cutoff rate specified,\nusing (# of points)/(space-time volume) =", k, ".\n")
	}
	if(k <= 0)
		stop("k must be greater than 0")
	num <- rpois(1, k*diff(X$xcoord)*diff(X$ycoord)*diff(X$tcoord))
	data <- data.frame(cbind(X$x, X$y, X$t))
	names(data) <- c("x", "y", "t")
	if(num > 0) {
		new.x <- runif(num, X$xcoord[1], X$xcoord[2]) 
		new.y <- runif(num, X$ycoord[1], X$ycoord[2])
		new.t <- runif(num, X$tcoord[1], X$tcoord[2])
		sim.pts <- data.frame(cbind(new.x,new.y,new.t))
		names(sim.pts) <- names(data) 
		sim.pts <- sim.pts[order(sim.pts$t), ]
		if(is.null(theta)) {
			fun2 <- function(z)
			{
				Y <- rbind(data, z)
				Y <- Y[which(Y$t <= z[3]),]
				Y <- stpp(Y$x, Y$y, Y$t, stwin(X$xcoord, X$ycoord, X$tcoord))
				tail(cifunction(Y),1)
			}
			lamb2 <- apply(sim.pts, 1, fun2)
		}
		if(!is.null(theta)) {
			fun2 <- function(z)
			{
				Y <- rbind(data, z)
				Y <- Y[which(Y$t <= z[3]),]
				Y <- stpp(Y$x, Y$y, Y$t, stwin(X$xcoord, X$ycoord, X$tcoord))
				tail(cifunction(Y, theta),1)
			}
			lamb2 <- apply(sim.pts, 1, fun2)
		}	
		prob <- (k-lamb2)/k
    	u <- runif(length(prob))
    	retain <- (u<=prob)
    	super <- data.frame(sim.pts[retain,])
    } else 
    	super <- data.frame()
	retain <- (lamb1 <= k)
	remove <- (lamb1 > k)
	keep1 <- data[retain, ]
	thin.data <- data[remove, ]
	thin.lamb <- lamb1[remove]
	prob <- k/thin.lamb
	u <- runif(length(prob))
	retain <- (u <= prob)
	keep2 <- thin.data[retain, ]
	deleted <- thin.data[!retain, ]
	allpts<-data.frame(rbind(keep1, keep2, super))
	y <- list(X = X, k = k, residuals = allpts, super = super, keep1 = keep1, keep2 = keep2, deleted = deleted)
	class(y) <- "superthin"
	return(y)
}
