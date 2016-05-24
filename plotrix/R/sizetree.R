sizetree<-function(x,left=0,top,right=1,lastcenter=NA,showval=TRUE,
 showcount=TRUE,stacklabels=TRUE,firstcall=TRUE,col=NULL,border=NA,
 toplab=NULL,base.cex=1,...) {

 dimx<-dim(x)
 colname<-names(x)[1]
 if(firstcall) {
  x<-x[do.call(order,x),]
  oldmar<-par("mar")
  par(mar=c(1,2,2,1))
  top<-sum(!is.na(x[,1]))
  if(top<dimx[1]) 
   cat(dimx[1]-top,"NA values dropped from first stack.\n")
  plot(0,xlim=c(0,dimx[2]),ylim=c(0,top),type="n", 
   axes=FALSE,xlab="",ylab="",...)
 }
 xfreq<-table(x[,1])
 lenxf<-length(xfreq)
 if(firstcall) {
  if(is.null(col)) {
   col<-list()
   for(index in 1:dimx[2])
    col[[index]]<-rainbow(length(table(x[,index])))
  }
  for(index in 1:dimx[2])
   if(is.null(names(col[[index]])))
    names(col[[index]])<-names(table(x[,index]))
 }
 if(lenxf) {
  if(is.list(col)) {
   barcol<-col[[1]]
   barcol<-barcol[names(col[[1]]) %in% names(xfreq)]
  }
  else barcol<-col[names(col) %in% names(xfreq)]
  labels<-names(xfreq)
  squeeze<-(right-left)/10
  for(bar in 1:lenxf) {
   if(length(xfreq[bar])) {
    if(!is.na(xfreq[bar])) {
     if(xfreq[bar] > 0) {
      rect(left+squeeze,top-xfreq[bar],right-squeeze,top,
       col=barcol[bar],border=border)
      labelheight<-strheight(labels[bar])
      cex<-ifelse((1.5*labelheight) > xfreq[bar], 
       base.cex*0.75*xfreq[bar]/labelheight,base.cex)
      if(showval) {
       textcol<-ifelse(colSums(col2rgb(unlist(barcol[bar])) *
        c(1.4,1.4,0.5)) < 350,"white","black")
       bartext<-ifelse(showcount,paste(labels[bar], 
        " (",xfreq[bar],")",sep=""),labels[bar])
       text((left+right)/2,top-xfreq[bar]/2, 
        bartext,cex=cex,col=textcol)
      }
      if(!is.na(lastcenter)) 
       segments(left+squeeze,top-xfreq[bar]/2,left-squeeze,
        lastcenter)
      xvalue<-ifelse(is.numeric(x[, 1]),as.numeric(labels[bar]),labels[bar])
      if(dimx[2] > 1) {
       newcol<-col
       newcol[[1]]<-NULL
       nextx<-subset(x,x[,1]==xvalue,2:dimx[2])
       sizetree(nextx,right,top,right+1,lastcenter=top-xfreq[bar]/2,
        showval=showval,firstcall=FALSE,col=newcol,border=border,base.cex=base.cex)
      }
     }
    }
   }
   top<-top-xfreq[bar]
  }
  if(stacklabels) mtext(colname,side=1,at=(left+right)/2,line=-0.4)
 }
 if(firstcall) {
  if(!is.null(toplab)) {
   par(xpd=TRUE)
   top<-sum(!is.na(x[,1]))
   text(0.5:(dimx[2]+0.5),1.01*top,toplab,adj=c(0.5,0))
   par(xpd=FALSE)
  }
  par(mar=oldmar)
 }
}
