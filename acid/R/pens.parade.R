pens.parade <-
function(x,bodies=TRUE,feet=0,...){
  x.wna<-x[!is.na(x)]
  n<- length(x.wna)
  x.sort<-sort(x.wna)
  xax<- seq(0,1,length=n)
  plot(xax,x.sort,pch=16,...)
  if(bodies==TRUE){
    for(i in 1:n){
      lines(c(xax[i],xax[i]),c(feet,x.sort[i]))
    }
  }
}
