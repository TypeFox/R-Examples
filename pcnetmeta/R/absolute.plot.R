absolute.plot <- function(nma.obj,alphabetic=TRUE,digits=2,width=5,height,network.name){
  if(!is.null(nma.obj$AbsoluteRisk)){
    ci<-nma.obj$AbsoluteRisk$Median_CI
    armparam<-"Absolute Risk"
  }
  if(!is.null(nma.obj$TrtEffect)){
    ci<-nma.obj$TrtEffect$Median_CI
    armparam<-"Treatment Effect"
  }
  if(!is.null(nma.obj$LogRate)){
    ci<-nma.obj$LogRate$Median_CI
    armparam<-"Log Rate"
  }
  ntrt<-dim(ci)[1]
  xx<-1:ntrt
  trtname<-rownames(ci)

  med<-low<-upp<-numeric(ntrt)
  for(i in 1:ntrt){
    str<-ci[i,1]
    split1<-strsplit(str,split=" \\(")
    med[i]<-as.numeric(split1[[1]][1])
    str2<-split1[[1]][2]
    split2<-strsplit(str2,split=", ")
    low[i]<-as.numeric(split2[[1]][1])
    upp[i]<-as.numeric(gsub("\\)","",split2[[1]][2]))
  }

  if(alphabetic){
    od<-order(trtname)
    trtname<-trtname[od]
    med<-med[od]
    low<-low[od]
    upp<-upp[od]
  }

  med.print<-format(round(med,digits=digits),nsmall=digits)
  low.print<-format(round(low,digits=digits),nsmall=digits)
  upp.print<-format(round(upp,digits=digits),nsmall=digits)
  cis<-paste(med.print," (",low.print,", ",upp.print,")",sep="")

  if(missing(height)) height<-ntrt-1
  if(missing(network.name)){
    filename<-paste("AbsolutePlot_",armparam,".pdf",sep="")
  }else{
    filename<-paste(network.name,"_AbsolutePlot_",armparam,".pdf",sep="")
  }
  pdf(filename,width=width,height=height)
  par(mfrow=c(1,1),mai=c(0.5,max(strwidth(trtname,"inches"))+0.1,0.1,max(strwidth(cis,"inches"))+0.1),mgp=c(1.5,0.5,0))
  plot(x=c(low,upp),y=rep(rev(1:ntrt),2),ylim=c(0.8,ntrt),frame.plot=FALSE,yaxt="n",xlab=armparam,ylab="",cex=0.1,col="white")
  points(med,rev(1:ntrt),pch=15)
  for(i in 1:ntrt){
    lines(x=c(low[i],upp[i]),y=c(ntrt+1-i,ntrt+1-i))
  }
  for(i in 1:ntrt){
    mtext(trtname[i],side=2,at=ntrt+1-i,las=1)
  }
  for(i in 1:ntrt){
    mtext(cis[i],side=4,at=ntrt+1-i,las=1)
  }
  par(mai=c(1.02,0.82,0.82,0.42),mgp=c(3,1,0))
  garbage<-dev.off()
  cat("The absolute plot is saved in your working directory:\n")
  cat(paste(getwd(),"\n",sep=""))
}
