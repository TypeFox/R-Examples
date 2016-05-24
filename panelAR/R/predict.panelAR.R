predict.panelAR <- function(object,newdata=NULL,se.fit = FALSE,
					conf.interval = FALSE, conf.level = 0.95, na.action=na.pass,...){
	beta <- object$coefficients
	tt <- terms(object)
	terms <- delete.response(tt)
	df <- object$df
	
	if (missing(newdata) || is.null(newdata)) {
        yhat <- object$fitted.values
        	 X <- model.matrix(terms,out$model)
    } else{
		mf <- model.frame(terms, newdata, na.action = na.action)
		# pull class descriptions of variables and see if they match the model frame
		.checkMFClasses(attr(terms, "dataClasses"), mf)
		X <- model.matrix(terms, mf)
		yhat <- X %*% beta
	}	

	if(se.fit || conf.interval){
		fit.out <- data.frame(fit=yhat)
		d <- diag(X%*%object$vcov%*%t(X))
		d[d<0] <- NA
		se <- sqrt(d)
		if(se.fit){
			fit.out$se <- se		
		}
		if(conf.interval){
			fit.out$lb <- yhat + se*qt((1-conf.level)/2,df=df)
        			fit.out$ub <- yhat + se*qt(1-(1-conf.level)/2,df=df)
		}
	} else{
		fit.out <- yhat
	}
	
out <- list(fit=fit.out,df=df)
out
}