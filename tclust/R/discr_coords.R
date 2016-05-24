discr_coords <-
function (x, equal.weights)
{
  if (is.null (x$par$x))
    stop ("dataset not included in tclust object - cannot calc discriminant coordinates.")

  if (missing (equal.weights))
    equal.weights = x$par$equal.weights

  x.no <- x$par$x[drop = FALSE, x$cluster != 0,]

  p <- ncol (x.no)
  n <- nrow (x.no)
  
  W <- matrix (0, ncol = p, nrow = p)
  cl.num = x$k # sum (x$usek)

  if (is.null (x$cov))		# tk-means items don't have cov info
  {
    if (is.null (x$par$x))
      stop ("Cannot access the data matrix. Set store.x = TRUE.")
	x$cov <- array (dim = c (p, p, cl.num))
	for (i in 1:cl.num)
		x$cov[,,i] <- cov (x$par$x[x$cluster == i,])
  }

  for (i in 1:cl.num)
  {
    if (equal.weights)
      W <- W + (n - 1) / cl.num * .getsubmatrix (x$cov, i)
    else
      W <- W + (x$size[i] - 1) * .getsubmatrix (x$cov, i)
  }

  wm <- eigen(W, symmetric = TRUE)
  wm$values[wm$values < 1e-06] <- 1e-06
  Tm <- t (wm$vectors %*% diag (sqrt (wm$values)))

  #Tm <- tdecomp.hh(W)            ##  
  Tinv <- solve (Tm)

  S <- (n - 1) * cov (x.no)
  B <- S - W
  Z <- t (Tinv) %*% B %*% Tinv
  dc <- eigen (Z, symmetric = TRUE)
  units <- Tinv %*% dc$vectors * sqrt (n - cl.num)
  proj <- x$par$x %*% units
  return (proj)
}
