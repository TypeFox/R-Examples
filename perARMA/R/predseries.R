predseries <-function(real,x,T,p,start,...)
{
 predseries_full<-function(real,x,T,p,start,missval,predcol,realcol){

pm<-permest(x,T, 0.05, missval,'series', pp=0)
pmean=pm$pmean
xd=pm$xd
n=length(xd)
pred<-predictperYW(xd,T,p,missval, start, predcol, realcol)
xp=pred$x

xpm<-matrix(0,start,1)
for(i in 1:start){
if (i%%T==0) {xpm[i]=xp[i]+pmean[T]
              } else {
             xpm[i]=xp[i]+pmean[i%%T]}
             }

plot (xpm[n:start], xlab="time",ylab= "values of the series", type="l", lwd=1, lty=1,col=predcol)
lines(real[n:start], type="l", lwd=1, lty=1,col=realcol)     
title(main="Prediction of the series vs. real data",
  sub = paste("observations from",length(x)," to", start))
legend("bottomright", c("real series","predicted series"), fill=c(realcol,predcol),ncol=2,title="legend")
}

L<-modifyList(list(missval=NaN,predcol="red",realcol="blue"), list(real=real,x = x, T=T, p=p,start=start,...))

 do.call(predseries_full,L)

 }
