Plot.PI <-
function(x,fore,Interval,start,end,frequency){
x=ts(x,start,end,frequency)

y1=ts(fore,start=end,frequency=frequency)
y2=ts(Interval,start=end,frequency=frequency)
ts.plot(x,y1,y2,lwd=c(1,rep(2,1+ncol(Interval))),col=c(1,4,rep(2,ncol(Interval))))
title(main="Time Plot and Prediction Intervals",sub="red = prediction quantiles; blue = point forecasts", col.sub=2)
d1=time(x); d2=d1[length(d1)]
abline(v=d2,col=3,lwd=2)
}
