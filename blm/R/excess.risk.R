crude.risk  <- function(formula, data, weights=NULL, na.action=na.omit){
	
risk <- function(start, stop, y, group, weights){
	lower.upper <- quantile(group, c(start, stop))
	index <- which(group>=lower.upper[1]&group<=lower.upper[2])
sum(y[index])/sum(weights[index])
}

covariate <- function(start, stop, group, weights){

	lower.upper <- quantile(group, c(start, stop))
	index <- which(group>=lower.upper[1]&group<=lower.upper[2])

sum(group[index]*weights[index])/sum(weights[index])
}

	data <- data[,all.vars(formula)]	
	if(is.null(weights)) weights <- rep(1, length(y))
	data$weights <- weights
	data <- na.action(data)
	
	y <- data[,all.vars(formula)[1]]
	group <- data[,all.vars(formula)[2]]
	
	Y <- mapply(risk, start=seq(0,.8,by=.01), stop=seq(.2,1,by=.01), MoreArgs=list(group=group, y=y, weights=data$weights))
	X <- mapply(covariate, start=seq(0,.8,by=.01), stop=seq(.2,1,by=.01), MoreArgs=list(group=group, weights=data$weights))

data.frame(risk = Y, x = X)
}


risk.exposure.plot <- function(object, scale=1,...){
	
	supplied.args <- list(...)
	
	if("ylab" %in% names(supplied.args))
		scatter.smooth(y=object$risk*scale, x=object$x, ...)
	else
		scatter.smooth(y=object$risk*scale, x=object$x, ylab="Average risk",...)
}