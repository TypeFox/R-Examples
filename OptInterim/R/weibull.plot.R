"weibull.plot" <-
function(param,x,l.type=1:3,l.col=c("blue","red"),...)
{  
  shape0 <- param[1]
  scale0 <- param[2]
  shape1 <- param[3]
  scale1 <- param[4]

  if(s(shape1,scale1,x)<=s(shape0,scale0,x))
  {
    stop("the alternative survial rate should be greater than the null survival rate")
  }
  
  ub <- 2*x
  a <- seq(0,2*x,length.out=1000)
  s0 <- s(shape0,scale0,a)
  s1 <- s(shape1,scale1,a)

  plot(a,s1,xlab="Time",ylab="Survival probability",main="Survival curves under the Weibull distribution", type="n",xlim=c(0,2*x),ylim=c(0,1),...)
  lines(a,s1,lty=l.type[1],col=l.col[1])
  lines(a,s0,lty=l.type[2],col=l.col[2])
  if(s(shape0,scale0,x/4)>0.4&s(shape0,scale0,3*x/4)>0.2)
    legend(x/4,0.2,legend=c("S1 (H1)","S0 (H0)"),lty=l.type[1:2],col=l.col,cex=0.8)
  else
  {  
    if(s(shape1,scale1,5*x/4)<0.7)
      legend(5*x/4,0.9,legend=c("S1 (H1)","S0 (H0)"),lty=l.type[1:2],col=l.col,cex=0.8)
    else
      legend(5*x/4,0.5,legend=c("S1 (H1)","S0 (H0)"),lty=l.type[1:2],col=l.col,cex=0.8)

  }  
  abline(v=x,lty=l.type[3])
  text(x,0.95*s(shape0,scale0,x),labels=round(s(shape0,scale0,x),3),cex=0.8,pos=2)
  text(1.1*x,1.05*s(shape1,scale1,x),labels=round(s(shape1,scale1,x),3),cex=0.8)
  text(1.07*x,0,labels=paste("x=",x,sep=" "),cex=0.8)
  
}
