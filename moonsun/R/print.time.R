`print.time` <-
function (x, ...)
 
{
	res = data.frame(time=as.vector(x));
	class(res$time)="time";
	rownames(res)=names(x);
	print(res);
	invisible(res);
}

