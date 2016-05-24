`as.gst` <-
function (x,jday=jd(),lambda=getOption("longitude"),...) 
{
	if (inherits(x,"gmt")) {

		S = jday - 2451545;
		t = S/36525;
		t0 = 6.697374558 + 2400.051336*t + 0.000025862*t^2;
		t0 = t0 %% 24;

		ut = x * 1.002737909;

		t0 = t0 + ut;
		t0 = t0 %% 24;

		t0[x == Inf] = Inf;
		t0[x == -Inf] = -Inf;

		class(t0) = c("gst","time");
		return(t0);
	}
	else
	if (inherits(x,"lst")) {

		lambda = lambda/15;
		res = x - lambda;
		res = res %% 24;

		res[x == Inf] = Inf;
		res[x == -Inf] = -Inf;


		class(res) = c("gst","time");
		return(res);
	}
	else
	if (inherits(x,"rst")) {
		res = x;
		res$rise=as.gst(x$rise,jday=jday,lambda=lambda,...);
		res$transit=as.gst(x$transit,jday=jday,lambda=lambda,...);
		res$set=as.gst(x$set,jday=jday,lambda=lambda,...);
		return(res);
	}
	else
	{
		res = as.vector(x);
		class(res) = c("gst","time");
		return(res);
	}		
}

