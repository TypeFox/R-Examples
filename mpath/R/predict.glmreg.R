predict.glmreg=function(object,newx,which=1:length(object$lambda), type=c("link","response","class", "coefficients","nonzero"), na.action=na.pass, ...){
 type=match.arg(type)
  if(missing(newx)){
    if(!match(type,c("coefficients","nonzero"),FALSE))stop("You need to supply a value for 'newx'")
     }
  else{
if(!is.null(object$terms)){
 mf <- model.frame(delete.response(object$terms), newx, na.action = na.action, xlev = object$xlevels)
 newx <- model.matrix(delete.response(object$terms), mf, contrasts = object$contrasts)
 newx <- newx[,-1] ### remove the intercept
}
}
  b0=as.matrix(object$b0)
  rownames(b0)="(Intercept)"
  nbeta=rbind2(b0,as.matrix(object$beta))[,which]
    vnames=dimnames(nbeta)[[1]]
    lambda=object$lambda[which]
  nlambda <- length(lambda)
  if(type=="coefficients")return(nbeta)
  if(type=="nonzero"){
  if(nlambda>1) return(nonzeroCoef(nbeta[-1,,drop=FALSE],bystep=TRUE))
  else return(which(abs(nbeta[-1]) > 0))
}
  famtype <- switch(object$family,
					"gaussian"=1,
					"binomial"=2,
					"poisson"=3,
					"negbin"=4)
  n <- dim(newx)[1]
  m <- dim(newx)[2]
  res <- .Fortran("pred",
				  n=as.integer(n),
				  m=as.integer(m),
				  nlambda=as.integer(nlambda),
				  x=as.double(newx),
				  b=as.double(object$beta[,which]),
				  a0=as.double(object$b0[which]),
                  family=as.integer(famtype),
				  eta = as.double(matrix(0,n,nlambda)),
				  mu = as.double(matrix(0,n,nlambda)),
				  package="mpath")
  eta <- matrix(res$eta, ncol=nlambda)
  mu <- matrix(res$mu, ncol=nlambda)
  colnames(eta) <- colnames(mu) <- colnames(object$beta[,which])
  #pihat <- exp(eta)/(1+exp(eta))
    if (object$family=="gaussian" | type=="link") return(eta)
    if (type=="response") return(mu)
	if (object$family=="binomial" & type=="class") return(eta>0)
  }
coef.glmreg <- function(object,which=1:length(object$lambda),...)
  {
  b0=matrix(object$b0, nrow=1)
  rownames(b0)="(Intercept)"
  nbeta=rbind2(b0,as.matrix(object$beta))
    return(nbeta[,which])
  }
