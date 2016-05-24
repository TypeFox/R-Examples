print.fads<-function(x,...) {
	UseMethod("print.fads")
}

print.fads.kfun<-function(x,...) {
	cat("Univariate second-order neighbourhood functions:\n")
	str(x)	
}

print.fads.k12fun<-function(x,...) {
	cat("Bivariate second-order neighbourhood functions:\n")
	str(x)	
}

print.fads.kpqfun<-function(x,...) {
	cat("Multivariate second-order neighbourhood functions :\n")
	cat("Interaction between each category p and each category q\n")
	str(x)
}

print.fads.kp.fun<-function(x,...) {
	cat("Multivariate second-order neighbourhood functions:\n")
	cat("Interaction between each category p and all the remaining categories.\n")
	str(x)
}

print.fads.kmfun<-function(x,...) {
	cat("Mark correlation functions:\n")
	str(x)
}

print.fads.ksfun<-function(x,...) {
	cat("Shimatani multivariate functions:\n")
	str(x)
}

print.fads.krfun<-function(x,...) {
	cat("Rao multivariate functions:\n")
	str(x)
}
