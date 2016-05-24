makeDendrite<-function(x) {
 dimx<-dim(x)
 if (is.null(dimx)) {
  dend.tab<-table(as.character(x),useNA="ifany")
  tablen<-length(dend.tab)
  dendrite<-vector("list",tablen)
  for(i in 1:tablen) dendrite[[i]]<-list(dend.tab[i],NULL)
 }
 else {
  dend.tab<- table(as.character(x[, 1]),useNA="ifany")
  tablen<-length(dend.tab)
  tabname<-names(dend.tab)
  dendrite<-vector("list",tablen)
  for(i in 1:tablen) {
   if(is.na(tabname[i]))
    nextx<-x[is.na(x[,1]),2:dimx[2]]
   else
    nextx<-x[as.character(x[,1])==tabname[i]&!is.na(x[,1]),2:dimx[2]]
   dendrite[[i]] <- list(dend.tab[i], makeDendrite(nextx))
  }
 }
 class(dendrite) <- "dendrite"
 return(dendrite)
}

sumDendrite<-function(x) {
 dsum<-0
 for(i in 1:length(x)) dsum<-dsum+x[[i]][[1]]
 return(dsum)
}

furc<-function(x,xpos,yrange,toplevel,maxx,cex=1,col) {
 xlen<-length(x)
 if(xlen) {
  yinc<-diff(yrange)/xlen
  ypos<-seq(yrange[1] + yinc/2, yrange[2] - yinc/2,length.out=xlen)
  if(xpos > maxx) xoffset<-rep(c(-0.17,0.17),length.out=xlen)
  else xoffset<-rep(0,xlen)
  if(!toplevel) {
   # only draw vertical segments if there is more than one category
   if(xlen > 1) segments(xpos-0.5,ypos[1],xpos-0.5,ypos[xlen])
   # horizontal segments
   segments(xpos-0.5,ypos,xpos+xoffset,ypos)
  }
  for(i in 1:xlen) {
   if(is.list(x[[i]][[2]])) {
    # draw the line to the next furcation
    segments(xpos,ypos[i],xpos+0.5,ypos[i])
    furc(x[[i]][[2]],xpos+1,c(ypos[i]-yinc/2, 
     ypos[i]+yinc/2),FALSE,maxx,cex=cex,col=col)
   }
   xname<-names(x[[i]][[1]])
   if(is.na(xname)) xname<-"NA"
   boxed.labels(xpos+xoffset[i],ypos[i],paste(xname,x[[i]][[1]]),
    cex=cex,bg=col[which(names(col)==xname)])
  }
 }
}

listDepth<-function(x) {
 if(is.list(x)) {
  maxdepth<-1
  for(lindex in 1:length(x)) {
   newdepth<-listDepth(x[[lindex]])+1
   if(newdepth > maxdepth) maxdepth<-newdepth
  }
 }
 else maxdepth<-0
 return(maxdepth)
}

plot.dendrite<-function(x,xlabels=NULL,main="",mar=c(1,0,3,0),cex=1,
 col="white",...) {

 if(class(x) != "dendrite") x<-makeDendrite(x)
 xnames<-unique(names(unlist(x)))
 xnames[is.na(xnames)]<-"NA"
 xnames<-sort(xnames)
 if(length(col) < length(xnames)) {
  col<-rep(col,length.out=length(xnames))
  names(col)<-as.character(xnames)
 }
 oldmar<-par("mar")
 par(mar=mar)
 xmax<-listDepth(x)/2
 ymax<-sumDendrite(x)
 plot(0,main=main,xlim=c(0.25,xmax),ylim=c(1,ymax), 
  xlab="",ylab="",type="n",axes=FALSE,...)
 par(xpd=TRUE)
 text(seq(0.5,xmax),par("usr")[3],xlabels)
 par(xpd=FALSE)
 furc(x,0.5,c(1,ymax),TRUE,maxx=xmax-1,cex=cex,col=col)
 par(mar=oldmar)
}

