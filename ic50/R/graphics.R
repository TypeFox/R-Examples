meanplot<-function(x,y,ic50,cl,cu,stddev,file,cpname){
  y_mean<-numeric(0)
  for(cn in 1:length(y[1,])) y_mean[cn]<-mean(y[,cn])
  plot(x,100*y_mean,xaxt="n",pch=".",
       xlab=expression(paste("concentration [",mu,"mol]",sep="")),ylab="viability [%]",
       xlim=c(floor(min(x)),ceiling(max(x))),ylim=c(0,max(100,max(100*(y_mean+stddev)))),
       main=paste("file ",strsplit(file,"/")[[1]][length(strsplit(file,"/")[[1]])]," (first), compound ",cpname,sep=""))
  range<-seq(-6,6,by=1)[floor(min(x))<=seq(-6,6,by=1) & seq(-6,6,by=1)<=ceiling(max(x))]
  axis(1,labels=c(as.character(10^range)),at=range)
  
  xgr<-linearCurve(x,y_mean)[[1]]
  ygr_mean<-linearCurve(x,y_mean)[[2]]
  lines(xgr,100*ygr_mean,lwd=2)

  for(cn in 1:length(y[1,])){
    lines(c(x[cn],x[cn]),100*c(y_mean[cn]-stddev[cn],y_mean[cn]+stddev[cn])) #Fehlerbalken
    lines(c(x[cn]-.025,x[cn]+.025),100*c(y_mean[cn]+stddev[cn],y_mean[cn]+stddev[cn]))
    lines(c(x[cn]-.025,x[cn]+.025),100*c(y_mean[cn]-stddev[cn],y_mean[cn]-stddev[cn]))
  }

  if(!is.na(ic50)) lines(c(ic50,ic50),c(-50,350))
  if(!is.na(cl) && !is.na(cu)){
    lines(c(cl,cl),c(-50,350),lty="dashed")
    lines(c(cu,cu),c(-50,350),lty="dashed")
  }
}



fittedplot<-function(x,y,ic50,cu,cl,file,cpname){
  y_mean<-numeric(0)
  for(cn in 1:length(y[1,])) y_mean[cn]<-mean(y[,cn])
  plot(x,100*y_mean,xaxt="n",
       xlab=expression(paste("concentration [",mu,"mol]",sep="")),ylab="viability [%]",
       xlim=c(floor(min(x)),ceiling(max(x))),ylim=c(0,max(100,max(100*y_mean))),
       main=paste("file ",strsplit(file,"/")[[1]][length(strsplit(file,"/")[[1]])]," (first), compound ",cpname,sep=""))
  range<-seq(-6,6,by=1)[floor(min(x))<=seq(-6,6,by=1) & seq(-6,6,by=1)<=ceiling(max(x))]
  axis(1,labels=c(as.character(10^range)),at=range)

  xgr<-linearCurve(x,y_mean)[[1]]
  ygr_mean<-linearCurve(x,y_mean)[[2]]
  lines(xgr,100*ygr_mean)

  if(!is.na(ic50)) lines(c(ic50,ic50),c(-50,350))
  if(!is.na(cl) && !is.na(cu)){
    lines(c(cl,cl),c(-50,350),lty="dashed")
    lines(c(cu,cu),c(-50,350),lty="dashed")
  }

  yexp<-expCurve(x,y_mean)[[2]]
  lines(xgr,100*yexp)
   
}



singleplot<-function(x,y,ic50,cu,cl,file,cpname){
  y_mean<-numeric(0)
  for(cn in 1:length(y[1,])) y_mean[cn]<-mean(y[,cn])
  
  plot(x,100*y_mean,xaxt="n",pch=".",
       xlab=expression(paste("concentration [",mu,"mol]",sep="")),ylab="viability [%]",
       xlim=c(floor(min(x)),ceiling(max(x))),ylim=c(0,max(100,max(100*y))),
       main=paste("file ",strsplit(file,"/")[[1]][length(strsplit(file,"/")[[1]])]," (first), compound ",cpname,sep=""))
  range<-seq(-6,6,by=1)[floor(min(x))<=seq(-6,6,by=1) & seq(-6,6,by=1)<=ceiling(max(x))]
  axis(1,labels=c(as.character(10^range)),at=range)

  lines(x,100*y_mean,lwd=2)
  
  if(!is.na(ic50)) lines(c(ic50,ic50),c(-50,350))
  if(!is.na(cl) && !is.na(cu)){
    lines(c(cl,cl),c(-50,350),lty="dashed")
    lines(c(cu,cu),c(-50,350),lty="dashed")
  }

  for(rep in 1:length(y[,1])){
    lines(x,100*y[rep,],lwd=2,lty="dotted")
  } 
}
