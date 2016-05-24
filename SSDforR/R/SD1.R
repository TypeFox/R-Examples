SD1 <-
function(behavior,phaseX,v1,ABxlab,ABylab, ABmain){
  
  maxy=which.max(behavior)
  max<-behavior[maxy]
  miny=which.min(behavior)
  min<-behavior[miny]
  
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  meanA=mean(A,na.rm=T)
  sdA=sd(A,na.rm=T)
  SDabove<-meanA+sdA
  SDbelow<-meanA-sdA
  #min=SDbelow-2
  #max=SDabove+2
  f1=SDabove >max
  f2=SDbelow <min
  
  if (f1==TRUE)
  {max=SDabove+1}
  
  
  if (f2==TRUE)
  {min=SDbelow-1}
  
  y<-na.omit(behavior)
  total=length(y)
  x=(1:total)
  
  end<-which(is.na(phaseX))
  np<-length(end)
  j=1
  while (j <= np){
    e<-end[j]
    
    
    
    y<-insert(y,NA,e)
    x<-insert(x,NA,e)
    j=j+1
  }
  
  graphics.off()
  layout(rbind(1,2), heights=c(6,1))
  plot(x,y,type="o",ylim=c(min,max),col="red",xlab=ABxlab,ylab=ABylab,main=ABmain,bty="l")
  #
  abline(h=meanA,col="green",lwd=3)
  abline(h=SDabove,col="black",lwd=3)
  abline(h=SDbelow,col="black",lwd=3)
  
  #par(mar=c(1, 1, 1, 1))
  #plot.new()
  #legend("center", c("behavior","+1sd","mean","-1sd"), col = c("red","black", "green","black"), lwd = 1,ncol=4,bty ="n")
  sdp<-c("SD",round(sdA,2))
  psdu<-c("+1SD",round(SDabove,2))
  pmean<-c("mean",round(meanA,2))
  psdb<-c("-1SD",round(SDbelow,2))
  
  
  tprint=c(sdp,psdu,pmean,psdb)
  print(tprint)
  ab<-NULL
  ab<<-recordPlot()
}
