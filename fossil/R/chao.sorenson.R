`chao.sorenson` <-
function(x,y)
{ n<-sum(x)
  m<-sum(y)
  if (length(x[y>0&x==2])==0) f2plus<-1
  else f2plus<-length(x[y>0&x==2])
  p1<-sum(x[y>0]/n)
  p2<-((m-1)/m)*(length(x[y>0&x==1])/(2*f2plus))
  p3<-sum(x[y==0]/n)
  u<-p1+p2*p3
  if (u>1) u<-1
  if (length(y[x>0&y==2])==0) fplus2<-1
  else fplus2<-length(y[x>0&y==2])
  q1<-sum(y[x>0]/m)
  q2<-((n-1)/n)*(length(y[x>0&y==1])/(2*fplus2))
  q3<-sum(y[x==0]/m)
  v<-q1+q2*q3
  if (v>1) v<-1
  if (u==0&v==0) c.s<-0
  else c.s<-(2*u*v)/(u+v)
  return(c.s)
  }

