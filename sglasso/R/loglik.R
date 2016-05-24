loglik <- function(object, N = 2){
	if(N <= 0) stop("N must be a positive integer")
	if(N != ceiling(N)) stop("N must be an integer value")
	p <- dim(object$S)[1]
	rho <- object$rho
	nv <- object$nv
	ne <- object$ne
	w <- object$w
	th_e <- object$theta[(nv + 1):(nv + ne), , drop = FALSE]
	l1_th_e <- apply(th_e, 2, function(x) {
					 id <- which(abs(x) > 0)
					 sum(w[id] * abs(x[id]))
					 }
					 )
	out <- Kh(object)
    out <- unlist(lapply(out, function(x) determinant(x, logarithm = TRUE)$modulus[1]))
	k <- 0.5 * N
	out <- k * (out - p + rho * l1_th_e)
	out
}
