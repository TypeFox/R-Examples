`gmt` <-
function (hour=NULL,minute=0,second=0,epoch=Sys.time(),length=1,by=1) 
{
	if (is.null(hour)) {
				daylt = as.POSIXlt(epoch,tz="GMT");
				hour = daylt$hour;
				minute = daylt$min;
				second = daylt$sec;
				}

	res = hour + minute/60 + second/3600;
	res = res %% 24;

	res = seq(res,length=length,by=by);

	class(res) = c("gmt","time");
	return(res);

}

