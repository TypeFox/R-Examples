##
##  b a r y c e n t r i c . R  Barycentric Lagrange Interpolation
##


barylag <- function(xi, yi, x) {
    stopifnot(is.vector(xi, mode="numeric"), is.vector(xi, mode="numeric"))
    if (!is.numeric(x))
        stop("Argument 'x' must be a numeric vector or matrix.")

	n <- length(xi); m <- length(x)

	# Check the input arguments
	if (length(yi) != n)
		stop("Node vectors xi an yi must be of same length.")
	if ( min(x) < min(xi) || max(x) > max(xi) )
		stop("Some interpolation points outside the nodes.")

	# Compute weights
	X  <- matrix(rep(xi, times=n), n, n)
	wi <- 1 / apply(X - t(X) + diag(1, n, n), 1, prod)

	# Distances between nodes and interpolation points
	Y <- outer(x, xi, "-")

	# Identify interpolation points that are nodes
	inds <- which(Y == 0, arr.ind=TRUE)
	Y[inds] <- NA

	# Compute the values of interpolated points
	W <- matrix(rep(wi, each=m), m, n) / Y
	y <- (W %*% yi) / apply(W, 1, sum)

	# Replace with values at corresponding nodes
	y[inds[,1]] <- yi[inds[,2]]

	# Return interpolation values as vector
	return(y[,])
}


barylag2d <- function(F, xn, yn, xf, yf) {
    M  <- nrow(F);     N <- ncol(F)
    Mf <- length(xf); Nf <- length(yf)

    # Compute weights
    X  <- matrix(xn, M, 1) %x% matrix(1, 1, M)  # Kronecker product
    Y  <- matrix(yn, N, 1) %x% matrix(1, 1, N)
    Wx <- t(1 / apply(X - t(X) + diag(1, M, M), 1, prod)) %x% matrix(1, Mf, 1)
    Wy <- t(1 / apply(Y - t(Y) + diag(1, N, N), 1, prod)) %x% matrix(1, Nf, 1)

    # Distances between nodes and interpolation points
    xdist <- matrix(xf, Mf, 1) %x% matrix(1,  1, M) -
             matrix(xn,  1, M) %x% matrix(1, Mf, 1)
    ydist <- matrix(yf, Nf, 1) %x% matrix(1,  1, N) -
             matrix(yn,  1, N) %x% matrix(1, Nf, 1)
             
    # Identify interpolation points that are nodes
    eps <- .Machine$double.eps
    xdist[xdist == 0] <- eps
    ydist[ydist == 0] <- eps  # Thanks to Greg von Winckel for this trick !

    Hx <- Wx / xdist
    Hy <- Wy / ydist

    Hx %*% F %*% t(Hy) /
        ( (matrix(apply(Hx, 1, sum), Mf, 1) %x% matrix(1, 1, Nf)) *
          (matrix(apply(Hy, 1, sum), 1, Nf) %x% matrix(1, Mf, 1)) )
}
