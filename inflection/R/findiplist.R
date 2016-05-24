findiplist <-
function(x,y,index)
{
  #Output for ESE, EDE methods as defined theoretically in:
  #Demetris T. Christopoulos(2012).'Developing methods for identifying the
  #inflection point of a convex/ concave curve'. arXiv:1206.5478v1 [math.NA]
  #Contact Emails: dchristop@econ.uoa.gr or dem.christop@gmail.com
  n=dim(x)[1]
  #For convex/concave data (upward sigmoid) give index=0
  #For concave/convex data (downward sigmoid) give index=1
  if(index==1){y=-y}
  x1=x[1]
  y1=y[1]
  x2=x[n]
  y2=y[n]
  flinp=function(x)
  {
    lin2(x1,y1,x2,y2,x)
  }
  if(n>=3)
  {
    LF=y-flinp(x)
    jf1=which.min(LF)
    xf1=x[jf1]
    jf2=which.max(LF)
    xf2=x[jf2]
    if(jf2<jf1)
    {xfx<-NaN}
    else
    {xfx<-.5*(xf1+xf2)}
  }
  else
  {jf1=NaN;jf2=NaN;xfx=NaN}
  if(n>=3)
  {
    vsl<-c()
    vsr<-c()
    for (i in(2:n)){
      A=findipl(x,y,i)
      vsl[i-1]<-A[3]
      vsr[i-1]<-A[4]
    }
    jl=which.min(vsl)+1
    xl=x[jl]
    jr=which.max(vsr)+1
    xr=x[jr]
    if((jl-jr>=2)==TRUE){xs<-.5*(xl+xr)}else{xs<-NaN;}
  }
  else
  {jl<-NaN;jr<-NaN;xs<-NaN;} 
  matrix(c(jr,jl,xs,jf1,jf2,xfx),nrow=2,ncol=3,byrow=TRUE)
}
