i.spline.x <-
function(xx,tt,i,k,delta=0.0001,Cs=F)
{
  nn=length(xx)
  a1=min(tt)
  b1=max(tt)
  xtt=seq(a1,b1,,by=delta)
  nxtt=length(xtt)
  ytt=double(nxtt)
  for(j in 1:nxtt)
  {
    
    ytt[j]=m.spline.x(xtt[j],tt,i,k)
  }
  ytt=cumsum(ytt)*delta
  if(Cs)
  {
    ytt=cumsum(ytt)*delta
  }
  
  ii=ceiling((xx-a1)/delta)
  ii[ii>nxtt]=nxtt
  ii[ii==0]=1
  res=ytt[ii]
  return(res)
}
