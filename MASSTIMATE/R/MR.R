MR <-
function(HC,FC,equation=c("raw","phylocor"),data=NULL) {
  if(equation=="raw") {
    log.MR<-1.78*log10(HC)+0.939*log10(FC)-0.215
    MR<-round(10^log.MR,1)
    ppe.err<-MR*0.24932
    upper.MR<-round(MR+ppe.err,1)
    lower.MR<-round(MR-ppe.err,1)
  }
  if(equation=="phylocor") {
    log.MR<-1.54*log10(HC)+1.195*log10(FC)-0.234
    MR<-round(10^log.MR,1)
    ppe.err<-MR*0.24624
    upper.MR<-round(MR+ppe.err,1)
    lower.MR<-round(MR-ppe.err,1)
  }    
  return(cbind(data,log.MR,MR,upper.MR,lower.MR))
}
