s2k <- function(y1,y2)
{
  if (length(y1)!=length(y2)) stop ("wrong number of elements")
  tablen<-length(y1)
  data<-c(y1,y2)
  x2a<-x2b<-0
  col1<-col2<-1
  p<-0
  z<-.C("x22k",data=as.integer(array(data)), tablen=as.integer(tablen),
        x2a=as.double(x2a), x2b=as.double(x2b), col1=as.integer(col1),
        col2=as.integer(col2), p=as.double(p),PACKAGE="gap")
  c1<-z$col1
  c2<-z$col2
  cat("\nthe maximum accumulated table below and above",c1,"\n")
  a<-sum(y1[1:c1])
  b<-sum(y1[-(1:c1)])
  c<-sum(y2[1:c1])
  d<-sum(y2[-(1:c1)])
  cat(a,b,a+b,"\n",c,d,c+d,"\n",a+c,b+d,a+b+c+d,"\n")
  cat("\nthe 1-to-other table with and without column",c2,"\n")
  a<-y1[c2]
  b<-sum(y1[-c2])
  c<-y2[c2]
  d<-sum(y2[-c2])
  cat(a,b,a+b,"\n",c,d,c+d,"\n",a+c,b+d,a+b+c+d,"\n")
  cat("\n")  

  list(x2a=z$x2a,x2b=z$x2b,col1=z$col1,col2=z$col2,p=z$p)
}
