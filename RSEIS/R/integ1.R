`integ1` <-
function (x,y,dm=-Inf,hm=+Inf)
{
  ####  integrate under the curve of a time series
  ####  return 2 numbers, one with the bottom triangle included
  ###  one without
  if(dm==-Inf)dm<-min(x)
  if(hm==+Inf)hm<-max(x)
  vyber<-x<=hm&x>=dm
  l<-length(x[vyber])
  v<-diff(x[vyber])
  z<-y[vyber][1:l-1]+y[vyber][2:l]
  o<-z*v/2
  osum<-sum(o)
  o1<-(y[x==min(x[vyber])]+y[x==max(x[vyber])])*(max(x[vyber])-min(x[vyber]))/2
  cista<-osum-o1
  return(c(osum,cista))
}

