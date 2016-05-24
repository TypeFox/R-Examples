summary.mine <- function(object, ...)
{
	result<-list()
	result$m<-object$m
	result$t<-object$t
	result$t.range<-object$t.range
	result$m.range<-object$m.range	

	class(result) <- "summary.mine"
	return(result)
}

print.mine <- function(x, ...)
{
	pp <- x$entropy.values
	colnames(pp) <- c("Time delay","Embedding dimension","Entropy")
	print(pp)

	invisible(x)
}

print.summary.mine <- function(x, ...)
{
	cat(paste("Embedding dimension: ",x$m, "[",paste(x$m.range,collapse=","),"]   \nTime delay: ",
	x$t, "[",paste(x$t.range,collapse=","),"]\n"));

	invisible(x)
}

