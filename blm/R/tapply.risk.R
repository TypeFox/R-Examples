tapply.risk <- function(formula,data,weights,values=NULL){
	
	y <- model.frame(formula,data)[,1]
	factor <- model.frame(formula,data)[,2]
	
	df <- data.frame(y=y,w=weights)
	df <- split(df,factor)
	
	risk <- sapply(df,function(x){sum(x$y*x$w)/sum(x$w)})
	
	if(is.null(values)){
		if(is.factor(values)) values <- levels(factor)
			else values <- sort(unique(factor))
			}

data.frame(risk=risk,x=values,n=as.numeric(table(factor)))
}