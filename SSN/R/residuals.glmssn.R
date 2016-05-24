residuals.glmssn <-
function(object, cross.validation=FALSE, ...)
{
	X <- object$sampinfo$X
	V <- object$estimates$V
	z <- object$sampinfo$z
	n <- object$sampinfo$sample.size - object$sampinfo$missing.sample.size
	p <- object$sampinfo$rankX
	betahat <- object$estimates$betahat
	data <- object$ssn.object@obspoints@SSNPoints[[1]]@point.data
	data1 <- data[object$sampinfo$ind.obs,]

	svd.out <- svd(V)
	V.5 <- svd.out$u %*% diag(1/sqrt(svd.out$d)) %*% t(svd.out$u)
	z.5 <- V.5 %*% z
	X.5 <- V.5 %*% X

	# raw residuals
	raw.resid <- as.vector(z - X %*% betahat)
	# Hat matrix, Montgomery and Peck, pg. 170
	Hat <- X.5 %*% invall(t(X.5) %*% X.5) %*% t(X.5)
	# standardized residuals, Montgomery and Peck, pg. 170
	e <- as.vector((diag(n) - Hat) %*% z.5)
	# studentized residuals Montgomery and Peck, pg. 172
	r <- as.vector(e/sqrt(diag((diag(n) - Hat))))
	# leverage - diagonal hat elements scaled by 2p/n, Montgomery and Peck, pg. 182
	hii <- diag(Hat)
	lever <- hii/(2*p/n)
	# Cook's D, Montgomery and Peck, pg. 182
	Di <- as.vector(r^2*hii/((1-hii)*p))
	outpt <- data.frame(obsval = z, pid = attr(z,"pid"), fitval = X %*% betahat, resid = raw.resid,
		resid.stand = e, resid.student = r, leverage = lever, CooksD = Di)
	names(outpt) <- c("obsval", "pid", "_fit_", "_resid_", "_resid.stand_",
		"_resid.student_", "_leverage_", "_CooksD_")
	if(cross.validation)
	{
		cv.out <- CrossValidationSSN(object)
		cv.resid <- as.vector(z - cv.out[,"cv.pred"])
		original.names <- names(outpt)
		outpt$resid.crossv <- cv.resid
		outpt$cv.pred <- cv.out$cv.pred
		outpt$cv.se <- cv.out$cv.se
		names(outpt) <- c(original.names, "_resid.crossv_", "_CrossValPred_", "_CrossValStdErr_")
	}
	data.tmp <- merge(data, outpt, by = "pid", all.x = TRUE)
  data <- data.tmp[order(data[,"pid"]),]
	object$ssn.object@obspoints@SSNPoints[[1]]@point.data <- data
	out <- object
	class(out) <- "influenceSSN"
	out
}

