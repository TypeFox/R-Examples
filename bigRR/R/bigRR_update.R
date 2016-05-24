bigRR_update <-
function(obj, Z, family = gaussian(link = identity), tol.err = 1e-6, tol.conv = 1e-8, GPU = FALSE)
{
	w <- as.numeric(obj$u^2/(1 - obj$leverage))
	w[w < tol.err] <- tol.err
	# Z.updated <- t(sqrt(w)*t(Z)) ## bug fixed 111201 -- xia
	bigRR(y = obj$y, X = obj$X, Z = Z, family = family, weight = w, tol.err = tol.err, tol.conv = tol.conv, GPU = GPU)
}