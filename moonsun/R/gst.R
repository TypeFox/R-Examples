`gst` <-
function (jday=jd(),hour=NULL,minute=0,second=0,epoch=Sys.time(),length=1,by=1) 
{
	S = jday - 2451545;
	t = S/36525;
	t0 = 6.697374558 + 2400.051336*t + 0.000025862*t^2;
	t0 = t0 %% 24;

	if (is.null(hour)) ut = gmt() else
	ut = hour + minute/60 + second/3600;

	ut = ut * 1.002737909;

	t0 = t0 + ut;
	t0 = t0 %% 24;

	t0 = seq(t0,length=length,by=by);

	class(t0) = c("gst","time");
	return(t0);

}

