`as.gmt` <-
function (x,jday=jd(),lambda=getOption("longitude"),...) 
{
	if (inherits(x,"gst")) {

		S = jday - 2451545;
		t = S/36525;
		t0 = 6.697374558 + 2400.051336*t + 0.000025862*t^2;
		t0 = t0 %% 24;

		res = x - t0;
		res = res %% 24;
		res = res * 0.9972695563;

		res[x == Inf] = Inf;
		res[x == -Inf] = -Inf;

		class(res) = c("gmt","time");
		return(res);
	}
	else
	if (inherits(x,"lt"))
       {
                d = as.numeric(as.POSIXct(format(Sys.time(),tz="GMT"))-as.POSIXct(format(Sys.time(),tz="")));
                res = x+d;
		    res = res %% 24;
                class(res)=c("gmt","time");
                return(res);
       }
	else
	if (inherits(x,"lst")) return(as.gmt(as.gst(x,lambda=lambda,...),jday=jday,...))
	else
	if (inherits(x,"rst")) {
		res = x;
		res$rise=as.gmt(x$rise,jday=jday,lambda=lambda,...);
		res$transit=as.gmt(x$transit,jday=jday,lambda=lambda,...);
		res$set=as.gmt(x$set,jday=jday,lambda=lambda,...);
		return(res);
	}
	else
	{
		res = as.vector(x);
		class(res) = c("gmt","time");
		return(res);
	}		
}

