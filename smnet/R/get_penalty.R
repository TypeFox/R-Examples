get_penalty<-function(ncol, cyclic = F){
  zvec<-rep(0, ncol)
  #   zvec[1]<-1;zvec[2]<- -1
  zvec[1]<--1;zvec[2]<- 2; zvec[3] <- -1
  circ<-circulant.spam(zvec)
  lastrw <- circ[nrow(circ),]
  circ   <- circ[-nrow(circ),]
  circ   <- circ[-nrow(circ),]
  out <- t(circ)%*%circ
  if(cyclic == T){out<-rbind(out, 1000*lastrw)}
  out
}
