XRchart <-
function(behavior,groupX,bandX,ABxlab,ABylab, ABmain){
  d2t=c(NA,1.128,1.693,2.059,2.326,2.534,2.704,2.847,2.970,3.078,3.173,3.259,3.336,3.407,3.472,3.532,3.588,3.640,3.689,3.735,3.778,3.819,3.858,3.895,3.931)  
  
  xm<-tapply(behavior,groupX,mean,na.rm=T)
  maxr<-tapply(behavior,groupX,max,na.rm=T)
  minr<-tapply(behavior,groupX,min,na.rm=T)
  p=xm
  r=maxr-minr
  minyw=which.min(xm)
  miny=xm[minyw]
  
  maxyw=which.max(xm)
  maxy=xm[maxyw]
  tn<-table(groupX)
  len<-length(tn)
  A<-seq(1:len)
  
  
  end<-which(is.na(behavior))
  lend=tn[1]
  d2=d2t[lend]
  e<-end[1]
  e=e-1
  e<-groupX[e]
  np<-length(end)
  A<-p[1:e]
  
  xr<-mean(r,na.rm=T)
  meanA=mean(xm,na.rm=T)
  sdA=(xr/d2)/sqrt(lend)
  
  Uband<-meanA+bandX*sdA
  Lband<-meanA-bandX*sdA
  x=(1:length(p))
  
  j=1
  while (j <= np){
    e<-end[j]
    e<-e-1
    e<-groupX[e]
    e<-e+j
    
    p<-insert(p,NA,e)
    x<-insert(x,NA,e)
    j=j+1
  }
  f1=Uband >maxy
  f2=Lband <miny
  
  if (f1==TRUE)
  {maxy=Uband+2}
  
  
  if (f2==TRUE)
  {miny=Lband-2}
  
  xmax<-length(p)
  
  graphics.off()
  plot.new()
  layout(rbind(1,2), heights=c(6,1))
  plot(x,p,ylim=c(miny,maxy),xlim=c(1,xmax),type="o",col="red",bty='L',main=ABmain,xlab=ABxlab,ylab=ABylab)
  
  abline(h=meanA,col="green")
  abline(h=Uband,col="blue")
  abline(h=Lband,col="orange")
  
  #par(mar=c(1, 1, 1, 1))
  #plot.new()
  #legend("center", c("behavior","Uband","mean","Lband"), col = c("red","blue", "green","orange"), lwd = 1,ncol=4,bty ="n")
  
  puband<-c("Uband=",round(Uband,3))
  pmean<-c("mean= ",round(meanA,3))
  plband<-c("Lband=",round(Lband,3))
  print(puband)
  print(pmean)
  print(plband)
  ab<-NULL
  
  ab<<-recordPlot()
}
