lambda0=function(x,y,weights=rep(1,N),exclude=NULL){
  if(length(exclude))x=x[,-exclude]
  N=length(y)

  ybar=weighted.mean(y,weights)
  yvar=weighted.mean((y-ybar)^2,weights)
  y=(y-ybar)/sqrt(yvar)
  weights=weights/N
  
  xbar=t(weights)%*%x
  xvar=t(weights)%*%(x^2)-xbar^2

  grad= abs(t(y*weights)%*%x)/sqrt(xvar)
  max(grad)
}
