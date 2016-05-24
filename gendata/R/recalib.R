recalib<-function(data,var,low,high){
  nmin<-low
  nmax<-high
  if(low>high){
    nmax<-low
    nmin<-high
  }
  cmin<-min(data[,var])
  cmax<-max(data[,var])
  
  data[,var]<-(nmax-nmin)/(cmax-cmin)*(data[,var]-cmin)+nmin
  return(data)
}
