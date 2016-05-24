select.glmnet.cr <-
function(fit, which="BIC") {
	x <- fit$x
	method <- fit$method
	predict.fit<-predict.glmnet.cr(fit,newx=x,method=method)
	if (c("BIC", "AIC")[charmatch(which, c("BIC","AIC"))] == "BIC") {which.min(predict.fit$BIC)}
	else if (c("BIC", "AIC")[charmatch(which, c("BIC","AIC"))] == "AIC") {which.min(predict.fit$AIC)}
}

