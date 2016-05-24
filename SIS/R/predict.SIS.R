predict.SIS <- function(object, newx, s = object$path.index, type = c("response","link","class"), ...){

if(class(object$fit)[1]=="cv.ncvreg"){ 
   family = object$fit$fit$family
   coefs = object$fit$fit$beta[-1,s]
   intercept = object$fit$fit$beta[1,s]
}else if(class(object$fit)[1]=="cv.glmnet"){ 
   family = switch(class(object$fit$glmnet.fit)[1], "elnet"="gaussian", "lognet"="binomial", "fishnet"="poisson", "coxnet"="cox")
   coefs = object$fit$glmnet.fit$beta[,s]
   intercept = object$fit$glmnet.fit$a0[s]
}else if(class(object$fit)[1]=="ncvreg"){ 
   family = object$fit$family
   coefs = object$fit$beta[-1,s]
   intercept = object$fit$beta[1,s]
}
else{
   family = switch(class(object$fit)[1], "elnet"="gaussian", "lognet"="binomial", "fishnet"="poisson", "coxnet"="cox")
   coefs = object$fit$beta[,s]
   intercept = object$fit$a0[s]
}

newx = newx[,object$ix]

if(length(coefs)==1 || length(newx)==1) eta = newx*coefs
else if(length(coefs)==length(s)) eta = newx%*%t(coefs)
else eta = newx%*%coefs

if(family!="cox") eta = t(intercept + t(eta))

if(family=="gaussian") pihat = eta
if(family=="binomial") pihat = exp(eta)/(1+exp(eta))
if(family%in%c("poisson","cox")) pihat = exp(eta)

if(match.arg(type)=="response") return(pihat)
if(match.arg(type)=="link") return(eta)
if(match.arg(type)=="class" && family=="binomial") return(matrix(as.integer(eta>0), ncol=length(s)))
if(match.arg(type)=="class" && family=="poisson") stop('Choose response for Poisson, or choose binomial for response')

}
