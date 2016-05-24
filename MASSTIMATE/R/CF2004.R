CF2004 <-
function(FC,data=NULL) {
  log.estimate<-2.738*log10(FC)-3.607
  CF2004<-round((10^log.estimate)*1000,2)
  return(cbind(data,CF2004))
}
