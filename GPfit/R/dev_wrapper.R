## Sets the parameter values equal and searches for a solution 
## with them equal
##
## May 8th, 2012

dev_wrapper <- function(beta,X,Y,nug_thres,...){
if (is.matrix(X) == FALSE){
	x = as.matrix(X)
}
d = ncol(X);
beta = rep(beta,d);
dev_val = GP_deviance(beta,X,Y,nug_thres,...);
return(dev_val)
}