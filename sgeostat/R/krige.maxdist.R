"krige.maxdist" <-
function (s, point.obj, at, var.mod.obj, maxdist) 
{
	    # Make a little distance function...
	    distance <- function(x1, y1, x2, y2) ((x1 - x2)^2 + (y1 - 
	    	    y2)^2)^0.5
	    if (!inherits(s, "point")) 
	    	    stop("prediction point s must be of class, \"point\".\n")
	    cat(paste("\nUsing points within", maxdist, "units of prediction points.\n", 
	    	    collapse = " "))
	    zhat <- NULL
	    sigma2hat <- NULL
	    cat("  Predicting")
	    for (looper in 1:(length(s$x))) {
	    	    cat(".")
	    	    x <- s$x[looper]
	    	    y <- s$y[looper]
	    	    #   calculate the distance between the prediction point and all points...
	    	    dst <- distance(point.obj$x, point.obj$y, x, y)
	    	    xvect <- point.obj$x[dst <= maxdist]
	    	    yvect <- point.obj$y[dst <= maxdist]
	    	    at2 <- at[dst <= maxdist]
	    	    dst <- dst[dst <= maxdist]
	    	    #    cat('length(xvect) = ',length(xvect),'\n')
	    	    if (length(xvect) == 0 | (!s$do[looper])) {
	    	    	    zhat <- c(zhat, NA)
	    	    	    sigma2hat <- c(sigma2hat, NA)
	    	    }
	    	    else {
	    	    	    #     Now construct the Big Gamma oh matrix...
	    	    	    distvect <- dist(cbind(xvect, yvect))
	    	    	    #     Distances come out as a vector, convert to full matrix...
	    	    	    n <- attr(distvect, "Size")
	    	    	    distmtrx <- matrix(0, n, n)
	    	    	    distmtrx[lower.tri(distmtrx)] <- distvect
	    	    	    distmtrx <- distmtrx + t(distmtrx)
	    	    	    GMatrix <- var.mod.obj$model(distmtrx, var.mod.obj$parameters)
	    	    	    nr <- if (is.null(nrow(GMatrix))) 
	    	    	    	    0
	    	    	    else nrow(GMatrix)
	    	    	    GMatrix <- cbind(GMatrix, rep(1, length = nr))
	    	    	    GMatrix <- rbind(GMatrix, c(rep(1, length = nr), 
	    	    	    	    0))
	    	    	    matrix.size <- length(xvect) + 1
	    	    	    #     Construct the little gamma oh vector...
	    	    	    gvector <- c(var.mod.obj$model(dst, var.mod.obj$parameters), 
	    	    	    	    1)
	    	    	    #     Solve the equations...
	    	    	    #      lambda.hat <- solve(qr(GMatrix),gvector)
	    	    	    #     Predict!
	    	    	    #      zhat <- c(zhat,sum(lambda.hat[1:(length(lambda.hat)-1)] * at2))
	    	    	    #      sigma2hat <- c(sigma2hat,sum(lambda.hat * gvector))
	    	    	    #
	    	    	    #     ... the above prediction mechanism would stop with an error
	    	    	    #         message, if the krige matrix is singular on some prediction 
	    	    	    #         point. The prediction on the next lines would not stop at all,
	    	    	    #         it only assigns NAs to such prediction points: 
	    	    	    #
	    	    	    Gqr <- qr(GMatrix)
	    	    	    if (Gqr$rank != 0) {
	    	    	    	    lambda.hat <- qr.coef(Gqr, gvector)
	    	    	    	    #     Predict!
	    	    	    	    zhat <- c(zhat, sum(lambda.hat[1:(length(lambda.hat) - 
	    	    	    	     1)] * at2))
	    	    	    	    sigma2hat <- c(sigma2hat, sum(lambda.hat * 
	    	    	    	     gvector))
	    	    	    }
	    	    	    else {
	    	    	    	    zhat <- c(zhat, NA)
	    	    	    	    sigma2hat <- c(sigma2hat, NA)
	    	    	    }
	    	    }
	    }
	    cat("\n")
	    #  return(point(s$x,s$y,list(zhat=zhat,sigma2hat=sigma2hat)))
	    s.o <- point(s)
	    s.o$zhat <- zhat
	    s.o$sigma2hat <- sigma2hat
	    return(s.o)
}
