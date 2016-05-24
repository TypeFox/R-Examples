`cusp.subspacerss` <-
function(predictors, dependents)
{
    X <- predictors
    Y <- dependents
	qx <- if(is.qr(X)) {X} else {qr(X)}
	qy <- if(is.qr(Y)) {Y} else {qr(Y)}
	dx <- qx$rank
	dy <- qy$rank
	#if(min(dx,dy)<2) browser()
	Qx <- qr.Q(qx)[,1:dx, drop=FALSE]
	Qy <- qr.Q(qy)[,1:dy, drop=FALSE]
	z <- svd(crossprod(Qx, Qy), nu=0)
	Ry <- qr.R(qy)[1:dy, 1:dy, drop=FALSE]
	rss <- (1-z$d^2) * colSums((t(Ry) %*% z$v)^2)
	list(rss = rss, cor=z$d)
}

