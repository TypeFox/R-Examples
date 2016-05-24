`anamorph` <-
function (x,k, plot=F) 
{
  theta<-hermite.coef(x,k)
  len<-length(x)
  gaus<-qnorm(seq(0,1,length=len+2)[-c(1,len+2)])
  x.g<-hermite(gaus,theta)
  x.g.s<-sort(x.g)
  #extrapolation
  x.g.low<-x.g[1:10]
  x.g.high<-x.g[(len-10):len]
  gaus.low<-gaus[1:10]
  gaus.high<-gaus[(len-10):len]
  low<-predict(lm(x.g.low~gaus.low),newdata=data.frame(gaus.low=-10))
  high<-predict(lm(x.g.high~gaus.high),newdata=data.frame(gaus.high=10))
  r<-range(x)
  if(plot){
  plot(c(-10,gaus,10),c(low,x.g,high), type="l", xlim=range(gaus)+c(-1,+1), ylim=r+diff(r)/4*c(-1,1))
  lines(c(-10,gaus,10),c(low,x.g.s,high), type="l", col=2)
  points(gaus,sort(x), col=2)
  }
  xtog<-approxfun(c(x.g.s,low,high),c(gaus,-10,10))
  gtox<-approxfun(c(gaus,-10,10),c(x.g.s,low,high))
  list(xtog=xtog,gtox=gtox)
}

