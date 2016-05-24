midks.plot <-
function(x.seq,y,dist,w.emp=NULL,...){
  n<-length(y)
  if(is.null(w.emp)) w.emp<-rep(1,n)
  w.emp<-w.emp/sum(w.emp)
  w.emp.cum<-cumsum(w.emp)
  cdens<-rep(NA,length(x.seq))
  for(i in 1:length(x.seq)){
    cdens[i]<-dist(x.seq[i],...)
  }
  plot(x.seq,cdens,type="l",ylim=c(0,1))
  x<-sort(y)
  for(i in 1:n){
    if(i==1){lines(c(0,x[i]),c(w.emp.cum[i],w.emp.cum[i]),col=3,lwd=2)
    }else{lines(c(x[i-1],x[i]),c(w.emp.cum[i],w.emp.cum[i]),col=3,lwd=2)} 
  }
}
