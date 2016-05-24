mBIC <-
function(b,data,N,L,h) {
  
  dif<-diff(b)
  
  if(N==0){
    mbic<-0
  } else if(any(dif<h)){
    mbic<- -100000
  } else {  
    
    mean.all<-mean(data[,1])
    ss.all<- sum((data[,1]-mean.all)^2)
    
    ### TERM 1
    segcal<-c()
    for(i in 1:(length(b)-1)) {
      
      if(dif[i]<=0) {
        dif[i]<-1
        segcal[i]<-0} 
      else{
        segmean<- sum(data[b[i]:(b[i+1]-1),1])/dif[i]     
        segcal[i]<-dif[i]*((segmean-mean.all)^2)
      }} 
    
    ss.bg<-sum(segcal)
    ss.wg<-ss.all-ss.bg
    
    term1<- 0.5*(L-N+1)*log(1+(ss.bg/ss.wg))
    
    ### TERM 2
    term2<-lgamma(0.5*(L-N+1))- lgamma(0.5*(L+1))
    
    ### TERM 3
    term3<-0.5*N*log(ss.all)
    
    ### TERM 4
    
    logdif<-log(dif)    
    term4<-0.5*sum(logdif)
    
    #### TERM 5
    term5<-(0.5-N)*log(L)
    
    mbic<- term1+term2+term3-term4+term5 } 
  
  return(mbic)
}
