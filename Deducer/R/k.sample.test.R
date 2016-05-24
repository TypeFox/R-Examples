k.sample.test<-function(formula,data=NULL,test=oneway.test,...){
	variables<-eval(formula[[2]],data,parent.frame())
	factor.var<-eval(formula[[3]],data,parent.frame())
	if(length(dim(variables))<1.5){
		variables<-d(variables)
		fn<-formula[[2]]
		names(variables)<-if(is.call(fn)) format(fn) else as.character(fn)
	}
	if(length(dim(factor.var))>1.5)
		factor.var <- factor.var[,1]
	

	x <- factor(factor.var)
	y <- variables
	tests <- list()
	tests[[1]] <- list(func = deparse(substitute(test)), level = levels(x))
	for (var.name in colnames(y)) {
		tmp<-na.omit(data.frame(y=y[,var.name],x=x))
		tests[[var.name]] <- test(y~x,tmp, ...)
		tests[[var.name]]$parameter = round(tests[[var.name]]$parameter,3)
	}
	result <- multi.test(tests)
	result
}