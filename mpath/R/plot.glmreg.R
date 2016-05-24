### from glmnet version 1.9-1 package
plot.glmreg=function(x, xvar=c("norm","lambda","dev"),label=FALSE,shade=TRUE, ...){
xvar=match.arg(xvar)
 plotCoef(x$beta,lambda=x$lambda,df=x$df,dev=x$dev.ratio,label=label,xvar=xvar,...)
}
