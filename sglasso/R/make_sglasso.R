make_sglasso <- function(object, call, algorithm, S, mask){
    nrho <- ifelse(is.null(object$nrho), 1, object$nrho)
	if(object$conv != 0){
		theta <- object$th[, 1:nrho, drop = FALSE]
		grd <- object$grd[, 1:nrho, drop = FALSE]
		rho <- object$rho[1:nrho]
	}	else {
		theta <- object$th
		grd <- object$grd
		rho <- object$rho
	}
    theta <- apply(theta, 2, function(x) zero(x, truncate = object$trnc))
    df <- apply(theta, 2, function(x) sum(abs(x) > 1.0e-13))
	w <- object$w
	flg <- as.logical(object$pnl_flg)
	nv <- object$nv
	ne <- object$ne
	obj <- list(call = call, nv = nv, ne = ne, theta = theta, w = w, flg = flg, df = df, rho = rho, grd = grd, 
				nstep = object$nstep, nrho = nrho, algorithm = algorithm, truncate = object$trnc, tol = object$tol, S = S,
				mask = mask, n = object$n, conv = object$conv)
	class(obj) <- "sglasso"
	obj
}
