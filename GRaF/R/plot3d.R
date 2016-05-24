plot3d <-
function (object, dims = c(1, 2), resolution = 20, CI = 0.95, prior = TRUE,
			CIs = TRUE, theta = 30, phi = 30, peak = FALSE, ...) {

#  facs <- which(unlist(lapply(object$obsx, is.factor)))
  facs <- object$facs
  if (any(dims %in% facs)) stop("3d plots are not implemented for factors")

  # set up matrix
  if (!peak) {
  	if(length(facs) > 0) {
		object$peak[(1:ncol(object$obsx))[-facs]] <- colMeans(object$obsx[, -facs, drop = FALSE])
		object$peak[, facs] <- sapply(data.frame(object$obsx[, facs, drop = FALSE]),
				function(x) names(sort(table(x), decreasing = TRUE))[1])
		for(i in facs) {
			object$peak[, i] <- factor(object$peak[, i], levels = levels(object$obsx[, i]))
		}
	} else {
		object$peak <- colMeans(object$obsx)
	}
  }


  pars <- object$peak
  k <- ncol(object$obsx)
  X <- as.data.frame(lapply(pars, rep, resolution ^ 2))
  seqx1 <- seq(min(object$obsx[, dims[1]]), max(object$obsx[, dims[1]]), length.out = resolution)
  seqx2 <- seq(min(object$obsx[, dims[2]]), max(object$obsx[, dims[2]]), length.out = resolution)
  expnd <- expand.grid(seqx1, seqx2)

  X[, dims[1]] <- expnd[, 1]
  X[, dims[2]] <- expnd[, 2]


#  scaledX <- X
#  notfacs <- (1:ncol(X))
#  if(length(facs) > 0) notfacs <- notfacs[-facs]
#  for (i in 1: length(notfacs)) {
#    scaledX[, notfacs[i]] <- (scaledX[, notfacs[i]] - object$scaling[1, i]) / object$scaling[2, i]
#  }

  # predict
  p <- predict(object, X, CI = CI)
  mn <- matrix(p[, 1], resolution)
  up <- matrix(p[, 3], resolution)
  low <- matrix(p[, 2], resolution)

  pri <- matrix(object$mnfun(X), resolution)

  # MAP
  persp(x = seqx1, y = seqx2, z = mn, zlim = c(0, 1), zlab = "probability of presence",
        xlab = colnames(object$x)[dims[1]], ylab = colnames(object$x)[dims[2]],
        col = 'light grey', border = 'dark grey', theta = theta, phi = phi,
        r = sqrt(20), ticktype = 'detailed', main = "posterior estimate", ...)
  if (prior) {
    persp(x = seqx1, y = seqx2, z = pri, zlim = c(0, 1), zlab = "probability of presence",
          xlab = colnames(object$x)[dims[1]], ylab = colnames(object$x)[dims[2]],
          col = 'light grey', border = 'dark grey', theta = theta, phi = phi,
          r = sqrt(20), ticktype = 'detailed', main = "prior", ...)
  }
  if (CIs) {
    persp(x = seqx1, y = seqx2, z = low, zlim = c(0, 1), zlab = "probability of presence",
          xlab = colnames(object$x)[dims[1]], ylab = colnames(object$x)[dims[2]],
          col = 'light grey', border = 'dark grey', theta = theta, phi = phi,
          r = sqrt(20), ticktype = 'detailed', main = "lower 95% CI", ...)
    persp(x = seqx1, y = seqx2, z = up, zlim = c(0, 1), zlab = "probability of presence",
          xlab = colnames(object$x)[dims[1]], ylab = colnames(object$x)[dims[2]],
          col = 'light grey', border = 'dark grey', theta = theta, phi = phi,
          r = sqrt(20), ticktype = 'detailed', main = "upper 95% CI", ...)
  }
}
