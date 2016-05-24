print.jd <-
function (x,...)
{
	res = data.frame(date=as.vector(x));
	class(res$date)="jd";
	rownames(res)=names(x);
	rrr <- data.frame(res, as.character(res))
	colnames(rrr) <- c("Date", "Julian Date") 
	print(rrr);
}

