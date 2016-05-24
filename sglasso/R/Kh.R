Kh <- function(object, rho){
	theta <- object$theta
	mask <- object$mask
	Dimnames <- dimnames(object$S)
	rho_obj <- object$rho
	if(!missing(rho)){
		if(!is.numeric(rho)){
			warning("rho is not numeric then it is set to the values used in sglasso")
			rho <- rho_obj
		}
		id <- rho_obj %in% rho
		if(!any(id)) stop("rho is not used in object")
		theta <- theta[, id, drop = FALSE]
		rho_obj <- rho_obj[id]
	}	
	out.Kh <- apply(theta, 2, function(th) theta2Kh(th, mask, Dimnames))
	names(out.Kh) <- paste("rho", rho_obj, sep = " = ")
	out.Kh
}