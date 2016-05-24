`as.lst` <-
function (x,lambda=getOption("longitude"),...) 
{

	if (inherits(x,"gst")) {

		if (is.null(lambda)) 
			{
			lambda = 0;
			warning("Your longitude is not set in the environment, assuming it is equal to 0");
			}

		res = x + lambda/15;
		res = res %% 24;

		res[x == Inf] = Inf;
		res[x == -Inf] = -Inf;

		names(res)=rownames(x);
		class(res)=c("lst","time");
		return(res);
	}
	else
	if (inherits(x,"gmt")) return(as.lst(as.gst(x,lambda=lambda,...),...))
	else
	if (inherits(x,"rst")) {
		res = x;
		res$rise=as.lst(x$rise,lambda=lambda,...);
		res$transit=as.lst(x$transit,lambda=lambda,...);
		res$set=as.lst(x$set,lambda=lambda,...);
		return(res);
	}
	else
	{
		res = as.vector(x);
		class(res)=c("lst","time");
		return(res);
	}		
}

