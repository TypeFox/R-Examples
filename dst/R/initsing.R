initsing<-function(card, names) { 
  # initialize results
  # card is the number of elements of the frame
  # names is a character vector of names
  x1<-diag(1,card)
  x2<-matrix(rep(1,times=card),nrow=1)
  f<-rbind(x1,x2)
  v<-matrix(c(rep(0,times=card),1),ncol=1)
  y<-bca(v,f,names)
  return(y)      
}