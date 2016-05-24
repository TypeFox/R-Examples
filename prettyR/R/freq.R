freq<-function(x,variable.labels=NULL,display.na=TRUE,decr.order=TRUE) {
 if(missing(x)) stop("A vector, dataframe or matrix must be supplied")
 xdim<-dim(x)
 if(is.null(xdim)) {
  if(is.null(variable.labels)) 
   variable.labels<-deparse(substitute(x))
  x<-list(x)
  nfreq<-1
 }
 else {
  nfreq<-xdim[2]
  if(is.matrix(x)) x<-as.data.frame(x)
  if(is.null(variable.labels)) variable.labels<-names(x)
  if(is.null(variable.labels)) variable.labels<-paste("V",1:xdim[2],sep="",collapse="")
 }
 freq.list<-rep(list(0),nfreq)
 for(i in 1:nfreq) {
  nna<-sum(is.na(x[[i]]))
  names(nna)<-"NA"
  freqs<-table(x[[i]])
  if(decr.order) freqorder<-order(freqs,decreasing=TRUE)
  else freqorder<-1:length(freqs)
  vl<-attr(x[[i]],"value.labels")
  if(!is.null(vl) & length(vl)) {
   vlabels<-attr(x[[i]],"value.labels")
   if(length(names(freqs)) < length(vlabels))
    vlabels<-vlabels[vlabels%in%names(freqs)]
   names(freqs)<-names(vlabels)
  }
  if(display.na) freqs<-c(freqs[freqorder],nna)
  else freqs<-freqs[freqorder]
  freq.list[[i]]<-freqs
 }
 names(freq.list)<-variable.labels
 class(freq.list)<-"freq"
 return(freq.list)
}
