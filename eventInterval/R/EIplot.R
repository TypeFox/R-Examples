EIplot<-function(event_times,main="",xlab="",ylab="",xaxticks=NA,xaxlabs=NA) {
 event_times<-sort(event_times)
 event_ints<-get_intervals(event_times)
 # drop the first event time
 event_times<-event_times[-1]
 plot(event_times,event_ints,type="b",main=main,xlab=xlab,ylab=ylab,xaxt="n")
 if(is.na(xaxticks[1])) axis(1)
 else axis(1,at=xaxticks,labels=xaxlabs)
 lines(supsmu(event_times,event_ints),lwd=2)
 rug(event_times)
}
