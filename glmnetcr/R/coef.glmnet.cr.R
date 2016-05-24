coef.glmnet.cr <-
function(object,s, ...) {
	list(a0=object$a0[s],beta=object$beta[,s])
}

