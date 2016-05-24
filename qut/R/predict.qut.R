predict.qut <- 
function(object, newx, mode='glm',offset=NULL,...){
	
	n=nrow(newx)
	
	#Check for warnings
	if(is.null(ncol(newx))) stop("newx should be a matrix with 2 or more columns")
	if(length(offset)!=n&!is.null(offset)) stop("length of offset is different to the number of observations")
	if(ncol(newx)!=length(object$beta[-1])) stop("newx should be a matrix with the same number of columns as covariates in the fitted model")
	
	if(mode=='lasso'){ #Predict with the lasso coefficients
		if(class(object$fit)[1]=='lars') type='fit' else type='response'
		newx=t(t(newx)/object$scale.factor)
		out=predict(object$fit,newx=newx,type=type,s=object$lambda,mode='lambda',offset=offset)
		if(class(object$fit)[1]=='lars') out=out$fit
	}
	else if(mode=='glm'){ #Predict with the glm fitted coefficients
		f=object$family()
		if(is.null(offset)) offset=0
		out=f$linkinv(object$betaglm[1]+newx%*%object$betaglm[-1]+offset)
	}
	
	return(out)
}