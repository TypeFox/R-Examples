ABiqr <-
function(behavior,phaseX,v1,ABxlab,ABylab, ABmain){
  
  maxy=which.max(behavior)
  max<-behavior[maxy]+2
  miny=which.min(behavior)
  min<-behavior[miny]-1
  
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  q=quantile(A,na.rm=T)
  
  
  p75<-q[4]
  p25<-q[2]
  medianA<-q[3]
  
  #min=p25-3
  #max=p25+3
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
  plot(x,y,ylim=c(min,max),type="o",col="red",xlab=ABxlab,ylab=ABylab,main=ABmain,bty="l")
  
  abline(h=medianA,col="green",lwd=3)
  abline(h=p75,col="blue",lwd=3)
  abline(h=p25,col="orange",lwd=3)
  
 # par(mar=c(1, 1, 1, 1))
  #plot.new()
  #legend("center", c("behavior","p75","median","p25"), col = c("red","blue", "green","orange"), lwd = 1,ncol=4,bty ="n")
  
  psdu<-c(round(p75,2))
  pmean<-c(round(medianA,2))
  psdb<-c(round(p25,2))
  iqr=q[4]-q[2]
  iqrp=c("IQR=",round(iqr,2))
  
  tprint=c(psdu,pmean,psdb)
  print(tprint)
  ab<-NULL
  ab<<-recordPlot()
}
