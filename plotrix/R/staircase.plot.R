staircase.plot<-function(heights,totals=NA,labels=NULL,halfwidth=0.3,main="",
 mar=NA,total.col="blue",inc.col=NA,bg.col=NA,direction="e",las=1,
 display.height=TRUE,stagger=FALSE,cex=par("cex"),prefix="",suffix="",...) {

 staircasePlot(heights=heights,totals=totals,labels=labels,halfwidth=halfwidth,
 main=main,mar=mar,stair.info=list(total.col=total.col,inc.col=inc.col,
 border=par("fg")),bg.col=bg.col,direction=direction,las=las,
 display.height=display.height,stagger=stagger,cex=cex,prefix=prefix,
 suffix=suffix,...)
}

staircasePlot<-function(heights,totals=NA,labels=NULL,halfwidth=0.3,main="",
 mar=NA,stair.info=list(total.col="blue",inc.col=NA,border=par("fg")),
 bg.col=NA,direction="e",las=1,display.height=TRUE,stagger=FALSE,cex=par("cex"),
 prefix="",suffix="",...) {

 if(is.matrix(heights) | is.data.frame(heights)) {
  dimheights<-dim(heights)
  if(dimheights[2] > 1) totals<-heights[,2]
  heights<-as.vector(heights)
 }
 if(!is.numeric(heights))
  stop("heights must be a numeric vector or matrix/data frame with numeric first column")
 nbars<-length(heights)
 if(length(prefix) < nbars) prefix<-rep(prefix,length.out=nbars)
 if(length(suffix) < nbars) suffix<-rep(suffix,length.out=nbars)
 # if no marker for increments (FALSE | 0) and totals (TRUE | non-zero)
 # consider the first and last values to be totals and all others increments 
 if(is.na(totals[1])) totals<-c(TRUE,rep(FALSE,nbars-2),TRUE)
 # coerce totals to a logical vector if it isn't
 if(!is.logical(totals)) totals<-totals != 0
 oldmar<-par("mar")
 if(!is.na(bg.col)) {
  oldbg<-par("bg")
  par(bg=bg.col)
 }
 maxht<-max(heights)
 # if there is a negative total height, make that the minimum, otherwise zero
 minht<-min(c(min(heights[totals]),0))
 currht<-heights[1]
 for(i in 2:nbars) {
  if(!totals[i]) {
   currht<-currht+heights[i]
   if(currht > maxht) maxht<-currht
  }
 }
 if(direction == "e" || direction == "w") {
  if(is.na(mar[1])) mar<-c(10,2,3,2)
  par(mar=mar,xaxs="i")
  plot(0,xlim=c(0.5,nbars+0.5),ylim=c(minht,maxht),type="n",axes=FALSE,
   xlab="",ylab="",...)
 }
 else {
  if(is.na(mar[1])) mar<-c(2,10,3,2)
  par(mar=mar,yaxs="i")
  plot(0,xlim=c(minht,maxht),ylim=c(0.5,nbars+0.5),type="n",axes=FALSE,
   xlab="",ylab="",...)
 }
 par(xpd=TRUE)
 bar.col<-rep(NA,nbars)
 if(length(stair.info$inc.col) < sum(!totals))
  stair.info$inc.col=rep(stair.info$inc.col,length.out=sum(!totals))
 bar.col[!totals]<-stair.info$inc.col
 bar.col[totals]<-stair.info$total.col
 label_offset<-ifelse(direction == "e" || direction == "w",
  strheight("M"),strwidth("M"))
 if(direction == "s" || direction == "w") {
  start<-nbars
  finish<-1
  dir<- -1
 }
 else {
  start<-1
  finish<-nbars
  dir<-1
 }
 barend<-0
 barpos<-start:finish
 for(bar in 1:nbars) {
  barstart<-ifelse(totals[bar],0,barend)
  barend<-barstart+heights[bar]
  if(direction == "e" || direction == "w") {
   rect(barpos[bar]-halfwidth,barstart,barpos[bar]+halfwidth,barend,
    col=bar.col[bar],border=stair.info$border)
   if(display.height)
    text(barpos[bar],ifelse(heights[bar]<0,barstart,barend)+label_offset,
     paste(prefix[bar],heights[bar],suffix[bar],sep=""),cex=cex)
   if(direction == "e" && bar != nbars)
    segments(barpos[bar]+halfwidth*dir,barend,barpos[bar]+dir-halfwidth*dir,
     barend,lty=3)
   if(direction == "w" && bar != nbars)
    segments(barpos[bar]-halfwidth*dir,barend,barpos[bar]+dir+halfwidth*dir,
     barend,lty=3)
   if(!is.null(labels)) {
    labelline<-0.2
    if(stagger) labelline<-c(labelline,1.2)
    mtext(labels,side=1,line=labelline,at=start:finish,adj=0.5,padj=1,
     las=las,cex=cex)
   }
  }
  else {
   rect(barstart,barpos[bar]-halfwidth,barend,barpos[bar]+halfwidth,
    col=bar.col[bar])
   if(display.height)
    text(ifelse(heights[bar]<0,barstart,barend)+label_offset,
     barpos[bar],paste(prefix[bar],heights[bar],suffix[bar],sep=""),adj=0,cex=cex)
   if(bar != nbars)
    segments(barend,barpos[bar]+halfwidth*dir,barend,
     barpos[bar]+dir-halfwidth*dir,lty=3)
   if(!is.null(labels)) {
    labelline<-0.5
    if(stagger) labelline<-c(labelline,1.5)
    mtext(labels,side=2,line=labelline,at=start:finish,adj=1,las=las,cex=cex)
   }
  }
 }
 if(nchar(main)) mtext(main,line=mar[3]/2,at=getFigCtr()[1],cex=1.5)
 par(xpd=FALSE,mar=oldmar)
 if(!is.na(bg.col)) par(bg=oldbg)
}
