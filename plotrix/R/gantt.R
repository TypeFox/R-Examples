get.gantt.info<-function(format="%Y/%m/%d") {
 cat("Enter the label, start and finish time for each task.\n")
 cat("Default format for time is year/month/day e.g. 2005/2/22\n")
 cat("Enter a blank label to end.\n")
 nextlabel<-"dummy"
 tasklabels<-NA
 taskstarts<-NA
 taskends<-NA
 priorities<-NA
 while(nchar(nextlabel)) {
  nextlabel<-readline("Task label - ")
  if(nchar(nextlabel)) {
   if(is.na(tasklabels[1])) tasklabels<-nextlabel
   else tasklabels<-c(tasklabels,nextlabel)
   nextstart<-as.POSIXct(strptime(readline("Task begins - "),format=format))
   if(is.na(taskstarts[1])) taskstarts<-nextstart
   else taskstarts<-c(taskstarts,nextstart)
   nextend<-nextstart-1
   while(nextend < nextstart) {
    nextend<-as.POSIXct(strptime(readline("Task ends - "),format=format))
    if(nextend < nextstart) cat("Task cannot end before it starts!\n")
    else {
     if(is.na(taskends[1])) taskends<-nextend
     else taskends<-c(taskends,nextend)
    }
   }
   nextpriority<-0
   while (nextpriority < 1 || nextpriority > 10) {
    nextpriority <- as.numeric(readline("Task priority (1-10) - "))
    if(is.na(nextpriority)) {
     cat("Task priority must be a number between 1 and 10!\n")
     next.priority<-0
    }
   }
   if(is.na(priorities[1])) priorities<-nextpriority
   else priorities<-c(priorities,nextpriority)
  }
 }
 return(list(labels=tasklabels,starts=taskstarts,ends=taskends,
  priorities=priorities))
}

gantt.chart<-function(x=NULL,format="%Y/%m/%d",xlim=NULL,taskcolors=NULL, 
 priority.legend=FALSE,vgridpos=NULL,vgridlab=NULL,vgrid.format="%Y/%m/%d",
 half.height=0.25,hgrid=FALSE,main="",xlab="",cylindrical=FALSE,label.cex=1,
 border.col=NA,priority.label="Priorities",priority.extremes=c("High","Low"),
 time.axis=3) {

 oldpar<-par("mai","omi","xpd","xaxs","yaxs")
 if(is.null(x)) x<-get.gantt.info(format=format)
 if(any(x$starts > x$ends)) stop("Can't have a start date after an end date")
 tasks<-unique(x$labels)
 ntasks<-length(tasks)
 if(is.null(x$priorities)) x$priorities <- rep(1, ntasks)
 if(is.null(dev.list())) plot.new()
 charheight<-strheight("M",units="inches")
 oldcex<-par(cex=label.cex)
 maxwidth<-max(strwidth(x$labels,units="inches")) + 0.3
 par(oldcex)
 if(is.null(xlim)) xlim = range(c(x$starts, x$ends))
 npriorities<-max(x$priorities)
 if(is.null(taskcolors)) taskcolors<-barcolors<-rainbow(npriorities)
 else barcolors<-taskcolors
 if(length(barcolors) < ntasks) barcolors <- barcolors[x$priorities]
 nlabels<-length(x$labels)
 bottom.margin<-ifelse(time.axis==3,0.8,1.3)
 if(length(barcolors) < nlabels) 
  barcolors<-barcolors[as.numeric(factor(x$labels))]
 par(mai=c(bottom.margin,maxwidth,charheight*5,0.1))
 par(omi=c(0.1,0.1,0.1,0.1),xaxs="i",yaxs="i")
 plot(range(x$starts),c(1,ntasks),xlim=xlim,ylim=c(0.5,ntasks+0.5),
  main="",xlab="",ylab="",axes=FALSE,type="n")
 box()
 if(nchar(main)) mtext(main,3,2,at=getFigCtr()[1])
 if(nchar(xlab)) mtext(xlab,1,1)
 if(is.na(vgrid.format)) {
  if(is.null(vgridlab)) vgridlab<-vgridpos
  axis(time.axis,at=vgridpos,labels=vgridlab,padj=0)
  tickpos<-vgridpos
 }
 else {
  if(is.null(vgridpos)) tickpos<-axis.POSIXct(3,xlim,format=vgrid.format)
  else tickpos<-vgridpos
  if(is.null(vgridlab) && !is.null(vgridpos)) 
   vgridlab<-format.POSIXct(vgridpos,vgrid.format)
  if(is.null(vgridlab)) axis.POSIXct(3,xlim,format=vgrid.format)
  else axis(time.axis,at=tickpos,labels=vgridlab,padj=0)
 }
 topdown<-seq(ntasks,1)
 axis(2,at=topdown,labels=tasks,las=2,cex.axis=label.cex)
 abline(v=tickpos,col="darkgray",lty=3)
 for(i in 1:ntasks) {
  if(cylindrical) 
   cylindrect(x$starts[x$labels == tasks[i]],topdown[i]-half.height,
    x$ends[x$labels == tasks[i]],topdown[i]+half.height,col=barcolors[i],
    gradient="y")
  else
   rect(x$starts[x$labels == tasks[i]],topdown[i]-half.height,
    x$ends[x$labels == tasks[i]],topdown[i]+half.height,
    col=barcolors[x$labels == tasks[i]],border=border.col)
 }
 if(hgrid) 
  abline(h=(topdown[1:(ntasks-1)]+topdown[2:ntasks])/2,col="darkgray",lty=3)
 if(priority.legend) {
  par(xpd=TRUE)
  plim<-par("usr")
  plot.width<-diff(plim[1:2])
  line.height<-strheight("W")
  color.legend(plim[1]+2*plot.width/5,
   plim[3]-(6.5-time.axis)*line.height,
   plim[1]+3*plot.width/5,
   plim[3]-(5.5-time.axis)*line.height,
   legend=c("High",rep("",length(taskcolors)-2),"Low"),
   rect.col=taskcolors)
  par(xpd=NA)
  text(plim[1]+plot.width/2,
   plim[3]-(7-time.axis)*line.height,
   priority.label,adj=c(0.5,1))
  par(xpd=FALSE)
 }
 par(oldpar)
 invisible(x)
}
