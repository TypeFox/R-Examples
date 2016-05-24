muvar <- function(n.loci=1,y1=c(0,1,1),y12=c(1,1,1,1,1,0,0,0,0),p1=0.99,p2=0.9)
{
  if(n.loci==1)
    .C("onelocus",y1=as.single(y1),p1=as.single(p1),PACKAGE="gap")
  else
    .C("twolocus",y12=as.single(y12),p1=as.single(p1),p2=as.single(p2),PACKAGE="gap")
}

