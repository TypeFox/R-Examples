"krige.all" <-
function (s, point.obj, at, var.mod.obj) 
{
	    # Make a little distance function...
	    distance <- function(x1, y1, x2, y2) ((x1 - x2)^2 + (y1 - 
	    	    y2)^2)^0.5
	    # Now construct the Big Gamma oh matrix...
	    cat("\nUsing all points.\n  Preparing the kriging system matrix...")
	    distvect <- dist(cbind(point.obj$x, point.obj$y))
	    # Distances come out as a vector, convert to full matrix...
	    n <- attr(distvect, "Size")
	    distmtrx <- matrix(0, n, n)
	    distmtrx[lower.tri(distmtrx)] <- distvect
	    distmtrx <- distmtrx + t(distmtrx)
	    GMatrix <- var.mod.obj$model(distmtrx, var.mod.obj$parameters)
	    #browser()
	    GMatrix <- cbind(GMatrix, rep(1, length = nrow(GMatrix)))
	    GMatrix <- rbind(GMatrix, c(rep(1, length = nrow(GMatrix)), 
	    	    0))
	    cat("\n  Inverting the matrix")
	    if (nrow(GMatrix) > 100) 
	    	    cat(" (and it's a big one)")
	    cat("...")
	    GMatrix.inv <- solve(qr(GMatrix))
	    cat("\n")
	    zhat <- NULL
	    sigma2hat <- NULL
	    cat("  Predicting")
	    for (looper in 1:(length(s$x))) {
	    	    cat(".")
	    	    x <- s$x[looper]
	    	    y <- s$y[looper]
	    	    if (!s$do[looper]) {
	    	    	    zhat <- c(zhat, NA)
	    	    	    sigma2hat <- c(sigma2hat, NA)
	    	    }
	    	    else {
	    	    	    #   calculate the distance between the prediction point and all points...
	    	    	    dst <- distance(point.obj$x, point.obj$y, x, y)
	    	    	    xvect <- point.obj$x
	    	    	    yvect <- point.obj$y
	    	    	    #   Construct the little gamma oh vector...
	    	    	    gvector <- c(var.mod.obj$model(dst, var.mod.obj$parameters), 
	    	    	    	    1)
	    	    	    #   Solve the equations...
	    	    	    lambda.hat <- GMatrix.inv %*% gvector
	    	    	    #   Predict!
	    	    	    zhat <- c(zhat, sum(lambda.hat[1:(length(lambda.hat) - 
	    	    	    	    1)] * at))
	    	    	    sigma2hat <- c(sigma2hat, sum(lambda.hat * gvector))
	    	    }
	    }
	    cat("\n")
	    #  return(point(s$x,s$y,list(zhat=zhat,sigma2hat=sigma2hat)))
	    s.o <- point(s)
	    s.o$zhat <- zhat
	    s.o$sigma2hat <- sigma2hat
	    return(s.o)
}
