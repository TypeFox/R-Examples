CrossValidationSSN <- function(object)
{
  z <- object$sampinfo$z
  X <- as.matrix(object$sampinfo$X)
  V <- object$estimates$V
  Vi <- object$estimates$Vi
  n <- object$sampinfo$obs.sample.size
  cdd.out <- matrix(-999.9, nrow = n, ncol = 3)
  cdd.out[,1] <- attributes(object$sampinfo$z)$pid

	for(i in 1:n) {
		Vi.i <- Vi[(1:n) != i,(1:n) != i] -
			matrix(Vi[(1:n) != i,i],ncol = 1) %*%
			matrix(Vi[i,(1:n) != i],nrow = 1)/Vi[i,i]
		c.i <- matrix(V[(1:n) != i,i],ncol = 1)
		xi <- matrix(X[i,], ncol = 1)
		X.i <- X[(1:n) != i,]
		z.i <- matrix(z[(1:n) != i], ncol = 1)
		xxi <- xi - t(X.i) %*% Vi.i %*% c.i
		covb.i <- solve(t(X.i) %*% Vi.i %*% X.i)
		si <- V[i,i]  - t(c.i) %*% Vi.i %*% c.i
		lam <- t(c.i + X.i %*% covb.i %*% xxi) %*% Vi.i

		cdd.out[i,2] <- lam %*% z.i
		cdd.out[i,3] <- sqrt(si + t(xxi) %*% covb.i %*% xxi)

	}
	cdd.out <- as.data.frame(cdd.out)
	names(cdd.out) <- c("pid","cv.pred","cv.se")
	cdd.out
}

