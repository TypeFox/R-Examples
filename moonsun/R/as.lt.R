`as.lt` <-
function (x,...) 
{
	if (inherits(x,"gmt"))
	{
		d = as.numeric(as.POSIXct(format(Sys.time(),tz="GMT"))-as.POSIXct(format(Sys.time(),tz="")));
		res = x-d;
		res = res %% 24;

		res[x == Inf] = Inf;
		res[x == -Inf] = -Inf;

		class(res)=c("lt","time");
		return(res);
	}
	else	if (inherits(x,"rst")) {
		res = x;
		res$rise=as.lt(x$rise,...);
		res$transit=as.lt(x$transit,...);
		res$set=as.lt(x$set,...);
		return(res);
	}
	else if (inherits(x,"lst")) 
		return(as.lt(as.gmt(x,...),...))
      else
        {
                res = as.vector(x);
                class(res) = c("lt","time");
                return(res);
        }  
}

