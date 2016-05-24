ppe <-
function(true, pred, abs = TRUE) {
  if(length(true)!=length(pred)) stop("Length of 'true' must equal 'pred'")
  if(abs==TRUE) {
    ppe<-abs(((true-pred)/pred)*100)
  }
  if(abs==FALSE) {
    ppe<-((true-pred)/pred)*100
  }
  mean.ppe<-mean(ppe,na.rm=TRUE)
  #calculate confidence intervals
    std<-sd(ppe,na.rm=TRUE)
    n<-length(ppe)
    ci.s.x<-sqrt((std^2)/n)
    t<-qt(0.975,df=n-1)
    ci.ppe<-c(lower95ci=mean.ppe-t*ci.s.x,upper95ci=mean.ppe+t*ci.s.x)
  range.ppe<-range(ppe,na.rm=TRUE)
  sd.ppe<-sd(ppe,na.rm=TRUE)
  result <- list(ppe.list=ppe,mean=mean.ppe,conf.int=ci.ppe,range=range.ppe,st.dev=sd.ppe)
  return(result)
}