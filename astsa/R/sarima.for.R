sarima.for <-
function(xdata,n.ahead,p,d,q,P=0,D=0,Q=0,S=-1,tol=sqrt(.Machine$double.eps),no.constant=FALSE){ 
  xname=deparse(substitute(xdata))
  xdata=as.ts(xdata) 
  n=length(xdata)
  constant=1:n
  xmean = rep(1,n);  if(no.constant==TRUE) xmean=NULL
  if (d==0 & D==0) {
    fitit=stats::arima(xdata, order=c(p,d,q), seasonal=list(order=c(P,D,Q), period=S),
            xreg=xmean,include.mean=FALSE, optim.control=list(reltol=tol));
    nureg=matrix(1,n.ahead,1)        
} else if (xor(d==1, D==1) & no.constant==FALSE) {
    fitit=stats::arima(xdata, order=c(p,d,q), seasonal=list(order=c(P,D,Q), period=S),
            xreg=constant,optim.control=list(reltol=tol));
    nureg=(n+1):(n+n.ahead)       
} else { fitit=stats::arima(xdata, order=c(p,d,q), seasonal=list(order=c(P,D,Q), period=S), 
            optim.control=list(reltol=tol));
          nureg=NULL   
}     
#--
 fore=stats::predict(fitit, n.ahead, newxreg=nureg)  
#-- graph:
  U = fore$pred + 2*fore$se
  L = fore$pred - 2*fore$se
   a=max(1,n-100)
  minx=min(xdata[a:n],L)
  maxx=max(xdata[a:n],U)
   t1=xy.coords(xdata, y = NULL)$x 
   if(length(t1)<101) strt=t1[1] else strt=t1[length(t1)-100]
   t2=xy.coords(fore$pred, y = NULL)$x 
   endd=t2[length(t2)]
   xllim=c(strt,endd)
  ts.plot(xdata,fore$pred,col=1:2, type="o", xlim=xllim, ylim=c(minx,maxx), ylab=xname) 
  lines(fore$pred, col="red", type="p")
  lines(U, col="blue", lty="dashed")
  lines(L, col="blue", lty="dashed") 
#
  return(fore)
}

