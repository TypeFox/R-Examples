axis.mult<-function(side=1,at=NULL,labels,mult=1,mult.label="",mult.line,
 mult.labelpos=NULL,...) {
 if(is.null(at)) at<-axTicks(side)
 if(missing(labels)) labels<-at/mult
 axis(side,at,labels,...)
 mult.label<-paste(mult.label," (* ",mult," )",sep="",collapse="")
 # multiplier position defaults to centered on the outside
 if(is.null(mult.labelpos)) mult.labelpos<-side
 edges<-par("usr")
 if(side %% 2) {
  # either top or bottom
  if(mult.labelpos %% 2) {
   adj<-0.5
   at<-(edges[1]+edges[2])/2
   if(missing(mult.line)) mult.line<-ifelse(mult.labelpos == side,3,0)
  }
  else {
   adj<-ifelse(mult.labelpos == 2,1,0) 
   at<-ifelse(mult.labelpos == 2,edges[1],edges[2])
   if(missing(mult.line)) mult.line<-1
  }
 }
 else {
  # either left or right
  if(mult.labelpos %% 2) {
   adj<-ifelse(mult.labelpos == 1,1,0) 
   at<-ifelse(mult.labelpos == 1,edges[3],edges[4])
   if(missing(mult.line)) mult.line<-1
  }
  else {
   adj<-0.5
   at<-(edges[3]+edges[4])/2
   if(missing(mult.line)) mult.line=ifelse(mult.labelpos == side,3,0)
  }
 }
 mtext(mult.label,side,mult.line,at=at,adj=adj,...)
}
