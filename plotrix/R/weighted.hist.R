get.breaks<-function(x,breaks) {
 # if a break computing function name is passed
 if(is.character(breaks)) 
  nbreaks<-do.call(paste("nclass",breaks,sep=".",collapse=""),list(x))
 # if breaks is numeric
 if(is.numeric(breaks)) {
  # if just the number of breaks is passed
  if(length(breaks) == 1) {
   nbreaks<-breaks
  }
  # otherwise assume that breaks specifies the breakpoints
  else return(breaks)
 }
 breakinc<-diff(range(x))/nbreaks
 breaks<-c(min(x),rep(breakinc,nbreaks))
 breaks<-cumsum(breaks)
 return(breaks)
}

weighted.hist<-function(x,w,breaks="Sturges",col=NULL,plot=TRUE,
 freq=TRUE,ylim=NA,ylab=NULL,xaxis=TRUE,...) {
 
 if(missing(x))
  stop("Usage: weighted.hist(x,...) vector of values x required")
 if(missing(w)) w<-rep(1,length(x))
 breaks<-get.breaks(x,breaks)
 width<-diff(breaks)
 diffx<-diff(range(x))
 equidist<-sum(width-width[1]) < diffx/1000
 nbreaks<-length(breaks)-1
 # make sure that the last break is greater than the maximum value
 lastbreak<-breaks[nbreaks+1]
 breaks[nbreaks+1]<-breaks[nbreaks+1]+diffx/1000
 if(diff(range(breaks)) < diffx)
   warning("Not all values will be included in the histogram")
 counts<-rep(0,nbreaks)
 for(bin in 1:nbreaks)
  counts[bin]<-sum(w[x >= breaks[bin] & x < breaks[bin+1]])
 density<-counts/sum(counts)
 if(freq) {
  if(is.null(ylab)) ylab<-"Frequency"
  heights<-counts
  if(!equidist) 
   warning("Areas will not relate to frequencies")
 }
 else {
  if(!equidist) {
   heights<-density*mean(width)/width
   heights<-heights/sum(heights)
  }
  else heights<-density
  if(is.null(ylab)) ylab<-"Density"
 }
 if(plot) {
  if(is.null(col)) col<-par("bg")
  if(is.na(ylim)) ylim<-c(0,1.1*max(heights,na.rm=TRUE))
  mids<-barplot(heights,width=width,col=col,space=0,ylim=ylim,ylab=ylab,...)
  tickpos<-c(mids-width/2,mids[length(mids)]+width[length(width)]/2)
  if(xaxis) axis(1,at=tickpos,labels=signif(c(breaks[1:nbreaks],lastbreak),3))
 }
 else mids<-breaks[-length(breaks)]+width/2
 invisible(list(breaks=breaks,counts=counts,density=density,
 mids=mids,xname=deparse(substitute(x)),equidist=equidist))
}
