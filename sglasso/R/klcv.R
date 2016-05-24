klcv <- function(object, X, scale = 1){
	if(missing(object))
	stop("object argument is missing")
	if(missing(X))
	stop("object argument is missing")
	if(!is.vector(scale))
	stop("scale is not a vector")
	if(length(scale) != 1)
	stop("scale is not a scalar")
	rho <- object$rho
	S <- as.matrix(object$S)
	p <- dim(S)[1]
	n <- dim(X)[1]
	K <- Kh(object)
	nrho <- length(K)
	K <- unlist(lapply(K, as.vector))
	K <- array(K, dim = c(p, p, nrho))
	storage.mode(n) <- "integer"
	storage.mode(p) <- "integer"
	storage.mode(X) <- "double"
	storage.mode(S) <- "double"
	storage.mode(nrho) <- "integer"
	storage.mode(K) <- "double"
	gdf <- double(nrho)
	info <- integer(1)
	out.fortran <- .Fortran("gdf_fun", n = n, p = p, X = X, S = S, nrho = nrho, Kh = K, gdf = gdf, info = info)
	if(out.fortran$info != 0) stop("error in dpotrf or dpotri subroutine")
	ll <- loglik(object, n)
	klcv.out <- - (ll - 0.5 * scale * out.fortran$gdf) / n
	pos <- which.min(klcv.out)
	min.klcv <- klcv.out[pos]
	rho.opt <- rho[pos]
	out <- list(klcv = klcv.out, rho = rho, loglik = ll, gdf = out.fortran$gdf, scale = scale, min.klcv = min.klcv, rho.opt = rho.opt, rhoid = pos)
	class(out) <- "klcv"
	out
}
