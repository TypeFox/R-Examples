plot.reco <- function(x, ...){
	#if(length(x)>1){stop("more than one reco model achieved. please provide just one")}
	x <- x[[1]]
	plot(R ~ Temp, data=x$model, ...)
	nd <- seq(min(x$model$Temp), max(x$model$Temp), length.out=100)
	lines(predict(x, newdata=data.frame(Temp=nd)) ~ nd)
}