Plot.Fore <-
function(x,fore,start,end,frequency)
{
x=ts(x,start,end,frequency)

f=ts(fore,start=end,frequency=frequency)

ts.plot(x,f,lwd=c(1,2),col=c(1,4))
title(main="Time Plot and Point Forecasts",sub="blue = point forecasts", col.sub=4)
d1=time(x); d2=d1[length(d1)]
abline(v=d2,col=3,lwd=2)
}
