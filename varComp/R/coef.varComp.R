coef.varComp <-
function(object, what=c('fixed','beta','random','varComp','var.ratio','tau'), ...)
{
# S3 method for getting coef's of class varComp object. 
# what:   'fixed' or 'beta': Return fixed effect parameters
#         'random' or 'varComp': Return REML estimates. 
#         'var.ratio' or 'tau': Return ratios of REML estimates to error variance. 
  what=match.arg(what)
  what=switch(what, fixed='beta', random='varComp', var.ratio='tau', what)

  if(what=='beta'){
	if(FALSE){
    call0=object$call
    call.args=names(call0)
    Y=if(is.null(object[['Y']])) eval(call0[['Y']], envir=object$frame) else object[['Y']]
	}
	if('fixef'%in%names(object)) return(object$fixef)
	if('Y'%in%names(object)) {
		Y = object$Y
	}else if ('model'%in%names(object)){
		Y = model.response(object$model)
	}else stop("response variable is not recored.")
    X=model.matrix(object, what='fixed')
    
	if(ncol(X)>0L){
		this.V=vcov(object, what='Y')
		this.Vbet=ginv(crossprod(X, solve(this.V,X)))
		this.bet=drop(this.Vbet%*%crossprod(X, solve(this.V, Y)))
		names(this.bet)=colnames(X)
	}else this.bet=numeric(0L)
    this.bet
  }else if (what=='varComp'){
    c(object$varComps, error=object$sigma2)
  }else if (what=='tau'){
    object$parms
  }
}

# fixef.varComp <-
# function(object, ...) coef(object, what='fixed', ...)
