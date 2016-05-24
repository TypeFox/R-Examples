# calculates the intervals between the times of the events
# most importantly, applies a correction for multiple events
# occurring within the same time increment
get_intervals<-function(event_times) {
 event_ints<-diff(event_times)
 # check for zero intervals, if so, spread each group out
 zerodiffs<-which(event_ints == 0)
 if(length(zerodiffs)) {
  start<-1
  while(start < length(event_times)) {
   if(event_times[start] == event_times[start+1]) {
    end<-start+1
    while(event_times[end] == event_times[end+1]) end<-end+1
    inc<-1/(1+end-start)
    seqend<-(end-start)/(1+end-start)
    event_times[start:end]<-event_times[start:end]+seq(0,seqend,by=inc)
    start<-end
   }
   else start<-start+1
  }
  event_ints<-diff(event_times)
 }
 return(event_ints)
}

# displays a histogram of the event intervals,
# a smoothed line of the histogram counts and
# a curve representing the best fit exponential distribution
EIdist<-function(event_times,nbreaks=10,main="",xlab="",ylab="",
 xaxticks=NA,xaxlabs=NA) {

 event_times<-sort(event_times)
 event_ints<-get_intervals(event_times)
 # drop the first event time
 event_times<-event_times[-1]
 event_density<-hist(event_ints,breaks=10,main=main,xlab=xlab,ylab=ylab)
 lines(supsmu(event_density$mids,event_density$counts),lwd=2)
 ei_fit<-fitdistr(event_ints,"exponential")
 lines(rescale(dexp(1:max(event_density$mids),rate=ei_fit$estimate),
  c(0,par("usr")[4])),lty=2)
}
