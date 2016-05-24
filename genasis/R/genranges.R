genranges<-function(x,delta=0) {
  range<-c((floor(((1+delta)*min(x,na.rm=TRUE)-delta*max(x,na.rm=TRUE))/(10^(floor(log10(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))-0)))*(10^(floor(log10(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))-0))),(ceiling(((1+delta)*max(x,na.rm=TRUE)-delta*min(x,na.rm=TRUE))/(10^(floor(log10(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))-0)))*(10^(floor(log10(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))-0))))
  return(range)
}